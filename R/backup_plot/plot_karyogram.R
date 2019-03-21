#' Karyogram-like chromosome overview
#'
#' Plot a karyogram-like chromosome overview with read counts and CNV-state
#' from a \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#'
#' @param model A \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#' @param file A PDF file where the plot will be saved.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
plotKaryogram <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL) {

    if (class(model)=='GRanges') {
        binned.data <- model
        model <- list()
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'), bhattacharyya=NA)
        plot.karyogram(model, both.strands=both.strands, file=file)
    } else if (class(model)=="aneuHMM") {
        plot.karyogram(model, both.strands=both.strands, file=file)
    } else if (class(model)=="aneuBiHMM") {
        plot.karyogram(model, both.strands=both.strands, plot.breakpoints=plot.breakpoints, file=file)
    }

}

# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
plot.karyogram <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL) {
    
    ## Check user input
    if (is.null(model$breakpoints) & plot.breakpoints) {
        warning("Cannot find breakpoint coordinates. Please run 'getBreakpoints' first.")
        plot.breakpoints <- FALSE
    }

    ## Convert to GRanges
  if (!is.null(model$bins$counts)) {
    bins <- model$bins
  } else if (!is.null(model$bincounts[[1]]$counts)) {
    bins <- model$bincounts[[1]]
    ind <- findOverlaps(bins, model$bins, select='first')
    bins$state <- model$bins$state[ind]
    bins$mstate <- model$bins$mstate[ind]
    bins$pstate <- model$bins$pstate[ind]
  }
    bins.split <- split(bins, seqnames(bins))
    bins.split <- bins.split[lengths(bins.split) > 0]

    ## Get some variables
    fs.x <- 13
    maxseqlength <- max(seqlengths(bins))
    tab <- table(bins$counts)
    tab <- tab[names(tab)!='0']
    if (both.strands) {
      custom.xlim <- get_rightxlim(c(bins$mcounts, bins$pcounts))
    } else {
      custom.xlim <- get_rightxlim(bins$counts)
    }

    # Quality info
    qualityInfo <- getQC(model)
    quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))

    ## Get breakpoint coordinates
    if (plot.breakpoints) {
        bp.coords <- model$breakpoints
        # Set to midpoint
        start(bp.coords) <- (start(bp.coords)+end(bp.coords))/2
        end(bp.coords) <- start(bp.coords)
    }

    ## Theme for plotting chromosomes
    empty_theme <- theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(size=fs.x),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

    ## Go through chromosomes and plot
    ggplts <- list()
    for (i1 in 1:length(bins.split)) {
      chrom <- names(bins.split)[i1]

        # Plot the read counts
        dfplot <- as.data.frame(bins.split[[i1]])
        # Transform coordinates to match p-arm on top
        dfplot$start <- (-dfplot$start + seqlengths(bins)[chrom])
        dfplot$end <- (-dfplot$end + seqlengths(bins)[chrom])
        # Set values too big for plotting to limit
        dfplot$counts[dfplot$counts>=custom.xlim] <- custom.xlim
        dfplot.points <- dfplot[dfplot$counts>=custom.xlim,]
        dfplot.points$counts <- rep(custom.xlim, nrow(dfplot.points))

        if (both.strands) {
            dfplot$mcounts <- - dfplot$mcounts    # negative minus counts
            dfplot$pcounts[dfplot$pcounts>=custom.xlim] <- custom.xlim
            dfplot$mcounts[dfplot$mcounts<=-custom.xlim] <- -custom.xlim
            dfplot.points.plus <- dfplot[dfplot$pcounts>=custom.xlim,]
            dfplot.points.plus$counts <- rep(custom.xlim, nrow(dfplot.points.plus))
            dfplot.points.minus <- dfplot[dfplot$mcounts<=-custom.xlim,]
            dfplot.points.minus$counts <- rep(-custom.xlim, nrow(dfplot.points.minus))
        }

        ## Read counts
        ggplt <- ggplot(dfplot, aes_string(x='start', y='counts'))    # data
        if (!is.null(bins.split[[i1]]$state)) {
            if (both.strands) {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts', col='pstate'), size=0.2)    # read count
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts', col='mstate'), size=0.2)    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts', col='pstate'), size=5, shape=21)    # outliers
                ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts', col='mstate'), size=5, shape=21)    # outliers
            } else {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts', col='state'), size=0.2)    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts', col='state'), size=2, shape=21)    # outliers
            }
          statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
            ggplt <- ggplt + scale_color_manual(values=stateColors(statelevels), drop=FALSE)    # do not drop levels if not present
        } else {
            if (both.strands) {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='pcounts'), size=0.2, col=strandColors('+'))    # read count
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mcounts'), size=0.2, col=strandColors('-'))    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')    # outliers
                ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='counts'), size=5, shape=21, col='gray20')    # outliers
            } else {
                ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts'), size=0.2, col='gray20')    # read count
                ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='counts'), size=2, shape=21, col='gray20')    # outliers
            }
        }
        ## Chromosome backbone
        if (both.strands) {
            ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=0, xmax=seqlengths(bins)[chrom], col='white', fill='gray20')    # chromosome backbone as simple rectangle
        } else {
            ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, xmax=seqlengths(bins)[chrom], col='white', fill='gray20')    # chromosome backbone as simple rectangle
        }
        if (plot.breakpoints) {
            df.bp <- as.data.frame(bp.coords[seqnames(bp.coords)==names(bins.split)[i1]])
            # Transform coordinates to match p-arm on top
            df.bp$start <- (-df.bp$start + seqlengths(bins)[chrom])
            df.bp$end <- (-df.bp$end + seqlengths(bins)[chrom])
            if (nrow(df.bp)>0) {
            statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
              suppressMessages( ggplt <- ggplt + scale_color_manual(values=c(breakpointColors(), stateColors(statelevels)), drop=FALSE) )
                ggplt <- ggplt + geom_segment(data=df.bp, aes_string(x='start', xend='start', color='type'), y=-custom.xlim, yend=-0.5*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'), alpha=0.5)
            }
        }
        ggplt <- ggplt + empty_theme    # no axes whatsoever
        ggplt <- ggplt + ylab(names(bins.split)[i1])    # chromosome names
        if (both.strands) {
            ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-custom.xlim,custom.xlim))    # set x- and y-limits
        } else {
            ggplt <- ggplt + coord_flip(xlim=c(0,maxseqlength), ylim=c(-0.6*custom.xlim,custom.xlim))    # set x- and y-limits
        }
        ggplts[[i1]] <- ggplt
        
    }
    names(ggplts) <- names(bins.split)
    
    ## Combine in one canvas
    fs.title <- 20
    nrows <- 2    # rows for plotting chromosomes
    nrows.text <- 2    # additional row for displaying ID and qualityInfo
    nrows.total <- nrows + nrows.text
    ncols <- ceiling(length(bins.split)/nrows)
    plotlist <- c(rep(list(NULL),ncols), ggplts, rep(list(NULL),ncols))
    cowplt <- plot_grid(plotlist=plotlist, nrow=nrows.total, rel_heights=c(2,21,21,2))
    cowplt <- cowplt + cowplot::draw_label(model$ID, x=0.5, y=0.99, vjust=1, hjust=0.5, size=fs.title)
    cowplt <- cowplt + cowplot::draw_label(quality.string, x=0.5, y=0.01, vjust=0, hjust=0.5, size=fs.x)
    if (!is.null(file)) {
        ggsave(file, cowplt, width=ncols*1.4, height=nrows*4.6)
    } else {
        return(cowplt)
    }

}
