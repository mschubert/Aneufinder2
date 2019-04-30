#plotProfile <- function(model) {
#    bins <- as.data.frame(model$bins)
#    if (!"midpoint" %in% colnames(bins))
#        bins$midpoint = (bins$start + bins$end) / 2
#
#    segs <- as.data.frame(model$segments)
#
#    ggplot() +
#        geom_point(data=bins, aes(x=midpoint, y=counts)) +
#        geom_segment(data=segs, aes(x=start, xend=end, y=mean.counts, yend=mean.counts), size=5) +
#        facet_grid(. ~ seqnames, scales="free_x") +
#        theme(panel.spacing = unit(0, "lines"),
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              panel.background = element_blank(),
#              axis.title.x = element_blank(),
#              axis.text.x = element_blank(),
#              axis.ticks.x = element_blank())
#}

plotProfile <- function(model) {
  # Profile (Top)
  bins                  <- model$bins
  plot.max.count        <- stats::quantile(bins$counts, 0.995)                                                         # Maximum count (xlim for every chromosome)
  if(length(plot.max.count) == 0 | is.na(plot.max.count) | is.infinite(plot.max.count)){                                          #  New code replacing function: 'get_rightxlim'
    plot.max.count         <- 1
  }
  bins                  <- transCoord(bins)
  dfplot                <- as.data.frame(bins)
  empty_theme           <- theme(axis.line=element_blank(),axis.ticks=element_blank(),
                                 axis.title.x=element_blank(),panel.background=element_blank(),
                                 panel.border=element_blank(),panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(),plot.background=element_blank())
  bprof                  <- ggplot(dfplot, aes_string(x='start.genome', y='counts'))                                    # ???    bpfl: bin profile
  bprof                  <- bprof + geom_jitter(aes_string(x='start.genome', y='counts'),
                                               position=position_jitter(width=0, height=0))                            # jitter is adding random variation change this line? replace by 'geom_point'?
  if (!is.null(model$segments$state)) {                                                                                # Are there cases where we don't have a copy number?
    dfplot.seg          <- as.data.frame(transCoord(model$segments))
    dfplot.seg$counts.CNV <- model$distributions[as.character(dfplot.seg$state),'mu'] 
    bprof               <- bprof + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',
                                   y='counts.CNV',xend='end.genome',yend='counts.CNV', col='state'), size=2)           # Line that is indicating the state
  }
  cum.seqlengths        <- cumsum(as.numeric(seqlengths(bins)))
  names(cum.seqlengths) <- seqlevels(bins)
  cum.seqlengths.0      <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- seqlevels(bins)
  df.chroms             <- data.frame(x=c(0,cum.seqlengths))
  bprof                 <- bprof + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2)     # Dashed lines between chromosomes
  bprof                 <- bprof + scale_color_manual(name="state", 
                                     values=stateColors(levels(dfplot.seg$state)), drop=FALSE)                         # do not drop levels if not present
  bprof                 <- bprof + empty_theme                                                                         # no axes whatsoever
  bprof                 <- bprof + coord_cartesian(ylim=c(0,plot.max.count))                                              # set x- and y-limits        # Should we not change the value of points higher than the limit to the limit????
  bprof                 <- bprof + scale_x_continuous(breaks=seqlengths(model$bins)/2+
                             cum.seqlengths.0[as.character(seqlevels(model$bins))], labels=seqlevels(model$bins))
  bprof                 <- bprof + ylab('read count') + ggtitle(model$ID) + theme(plot.title = element_text(hjust = 0.5))

  # Histogram (Bottom)
  counts                <- bins$counts
  bin.width             <- round(plot.max.count / 50)
  if(bin.width < 1){
    bin.width           <- 1
  }
  bhist                 <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'),
                                  binwidth=bin.width, color='black', fill='white') + coord_cartesian(xlim=c(0,plot.max.count)) +
                                  theme_bw() + xlab("read count")
  if (!is.null(model$bins$state) & !is.null(model$weights)) {
    weights             <- model$weights
    c.state.labels      <- as.character(levels(model$bins$state))                                                      # Add fits to the histogram
    numstates           <- length(weights)
    x                   <- 0:max(counts)
    distributions       <- data.frame(x)
    for (istate in 1:nrow(model$distributions)) {
      if (model$distributions[istate,'type']=='delta') {
        distributions[[length(distributions)+1]] <- c(weights[istate],rep(0,length(x)-1))                              # zero-inflation
      } else if (model$distributions[istate,'type']=='dgeom') {
        distributions[[length(distributions)+1]] <- weights[istate] * stats::dgeom(                                    # geometric
                                                      x, model$distributions[istate,'prob'])
      } else if (model$distributions[istate,'type']=='dnbinom') {
        distributions[[length(distributions)+1]] <- weights[istate] * stats::dnbinom(x,                                # negative binomials
                                               model$distributions[istate,'size'], model$distributions[istate,'prob'])
      } else if (model$distributions[istate,'type']=='dpois') {
        distributions[[length(distributions)+1]] <- weights[istate] * stats::dpois(                                    # poissons
                                                      x, model$distributions[istate,'prob'])
      } else if (model$distributions[istate,'type']=='dbinom') {
        s               <- model$distributions[istate,'size']                                                          # binomials
        p               <- model$distributions[istate,'prob']
        distributions[[length(distributions)+1]] <- weights[istate] * stats::dbinom(x, round(s), p)
      }
    }
    distributions       <- as.data.frame(distributions)
    names(distributions) <- c("x",c.state.labels)
    distributions$total <- apply(distributions[,-1], 1, sum)                                                           # Total
    distributions       <- stats::reshape(distributions, direction="long", varying=1+1:(numstates+1),                  # Reshape the data.frame for plotting with ggplot
                                          v.names="density", timevar="state", times=c(c.state.labels,"total"))
    bhist               <- bhist + geom_line(aes_string(x='x', y='density', group='state', col='state'),               # Plot the distributions
                                             data=distributions)
    lmeans              <- round(model$distributions[,'mu'], 2)                                                        # Make legend and colors correct
    lvars               <- round(model$distributions[,'variance'], 2)
    lweights            <- round(model$weights, 2)
    legend              <- paste0(c.state.labels, ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
    legend              <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=',round(var(counts),2)))
    bhist               <- bhist + scale_color_manual(breaks=c(c.state.labels, 'total'),values=stateColors(c(c.state.labels,'total')), labels=legend)                                 
  }
  cowplt                <- cowplot::plot_grid(bprof, bhist, nrow=2, rel_heights=c(1.2,1))
  return(cowplt)
}
