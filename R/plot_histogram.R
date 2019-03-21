plotBivariateHistograms <- function(bihmm) {

    ## Stack the two strands (bins)
    binned.data <- bihmm$bins
    
    binned.data.minus <- binned.data
    strand(binned.data.minus) <- '-'
    mcols(binned.data.minus) <- NULL
    binned.data.minus$state <- binned.data$mstate
    binned.data.minus$copy.number <- binned.data$mcopy.number
    
    binned.data.plus <- binned.data
    strand(binned.data.plus) <- '+'
    mcols(binned.data.plus) <- NULL
    binned.data.plus$state <- binned.data$pstate
    binned.data.plus$copy.number <- binned.data$pcopy.number
    
    binned.data.stacked <- c(binned.data.minus, binned.data.plus)
    
    ## Attributes
    mask.attributes <- c('complexity', 'spikiness', 'entropy')
    attributes(binned.data.stacked)[mask.attributes] <- attributes(binned.data)[mask.attributes]
    
    ## Stack the two strands (bincounts)
    bincounts <- bihmm$bincounts[[1]]
    
    bincounts.minus <- bincounts
    strand(bincounts.minus) <- '-'
    mcols(bincounts.minus) <- NULL
    bincounts.minus$counts <- bincounts$mcounts
    
    bincounts.plus <- bincounts
    strand(bincounts.plus) <- '+'
    mcols(bincounts.plus) <- NULL
    bincounts.plus$counts <- bincounts$pcounts
    
    bincounts.stacked <- c(bincounts.minus, bincounts.plus)

    ## Make fake uni.hmm and plot
    strand <- 'minus'
    uni.hmm <- list()
    uni.hmm$ID <- bihmm$ID
    uni.hmm$bins <- binned.data.stacked
    uni.hmm$bincounts <- GRangesList(bincounts.stacked)
    uni.hmm$segments <- bihmm$segments
    uni.hmm$weights <- bihmm$univariateParams$weights
    uni.hmm$distributions <- bihmm$distributions[[strand]]
    uni.hmm$qualityInfo <- bihmm$qualityInfo
    class(uni.hmm) <- "aneuHMM"
    ggplts <- plotHistogram(uni.hmm)

    return(ggplts)

}

#' Plot a histogram of binned read counts with fitted mixture distribution
#'
#' Plot a histogram of binned read counts from with fitted mixture
#' distributions from a \code{\link{aneuHMM}} object.
#'
#' @param model A \code{\link{aneuHMM}} object.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom stats dgeom dnbinom dbinom dpois reshape
plotHistogram <- function(model) {

    model <- suppressMessages( loadFromFiles(model, check.class=c('GRanges', 'GRangesList', "aneuHMM"))[[1]] )
    if (class(model) == 'GRanges') {
        bins <- model
        bincounts <- model
    } else if (is(model, "GRangesList")) {
        bins <- model[[1]]
        bincounts <- model[[1]]
    } else if (class(model) == "aneuHMM") {
        bins <- model$bins
        if (!is.null(model$bins$counts)) {
            bincounts <- model$bins
        } else if (!is.null(model$bincounts[[1]]$counts)) {
            bincounts <- model$bincounts[[1]]
        }
    }
    counts <- bincounts$counts

      states <- bins$state
        if (!is.null(model$weights)) {
              weights <- model$weights
        }

      # Find the x limits
      breaks <- max(counts)
      if (max(counts)==0) { breaks <- 1 }
      rightxlim <- get_rightxlim(counts)

      # Quality info
      qualityInfo <- getQC(model)
      quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))

      # Plot the histogram
      ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count") + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
      if (is.null(model$weights)) {
            return(ggplt)
      }

      if (!is.null(model$bins$state)) {
          ### Add fits to the histogram
          c.state.labels <- as.character(levels(model$bins$state))
          numstates <- length(weights)
          x <- 0:max(counts)
          distributions <- data.frame(x)
      
          for (istate in 1:nrow(model$distributions)) {
                if (model$distributions[istate,'type']=='delta') {
                      # zero-inflation
                      distributions[[length(distributions)+1]] <- c(weights[istate],rep(0,length(x)-1))
                } else if (model$distributions[istate,'type']=='dgeom') {
                      # geometric
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dgeom(x, model$distributions[istate,'prob'])
                } else if (model$distributions[istate,'type']=='dnbinom') {
                      # negative binomials
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dnbinom(x, model$distributions[istate,'size'], model$distributions[istate,'prob'])
                } else if (model$distributions[istate,'type']=='dpois') {
                      # poissons
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dpois(x, model$distributions[istate,'prob'])
                } else if (model$distributions[istate,'type']=='dbinom') {
                      # binomials
                      s <- model$distributions[istate,'size']
                      p <- model$distributions[istate,'prob']
                      distributions[[length(distributions)+1]] <- weights[istate] * stats::dbinom(x, round(s), p)
                }
          }
          distributions <- as.data.frame(distributions)
          names(distributions) <- c("x",c.state.labels)
          # Total
          distributions$total <- apply(distributions[-1], 1, sum)
      
          # Reshape the data.frame for plotting with ggplot
          distributions <- stats::reshape(distributions, direction="long", varying=1+1:(numstates+1), v.names="density", timevar="state", times=c(c.state.labels,"total"))
          ### Plot the distributions
            ggplt <- ggplt + geom_line(aes_string(x='x', y='density', group='state', col='state'), data=distributions)
          
          # Make legend and colors correct
          lmeans <- round(model$distributions[,'mu'], 2)
          lvars <- round(model$distributions[,'variance'], 2)
          lweights <- round(model$weights, 2)
          legend <- paste0(c.state.labels, ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
          legend <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=', round(var(counts),2), ', weight(data)=1'))
          ggplt <- ggplt + scale_color_manual(breaks=c(c.state.labels, 'total'), values=stateColors(c(c.state.labels,'total')), labels=legend)
      }
  
      return(ggplt)

}
