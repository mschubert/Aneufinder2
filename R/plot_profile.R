#' Read count and CNV profile
#'
#' Plot a profile with read counts and CNV-state from a \code{\link{aneuHMM}}
#' object or \code{\link{binned.data}}.
#'
#' @param model A \code{\link{aneuHMM}} object or \code{\link{binned.data}}.
#' @param file A PDF file where the plot will be saved.
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @param both.strands If \code{TRUE}, strands will be plotted separately.
#' @param normalize.counts An character giving the copy number state to which
#'   to normalize the counts, e.g. '1-somy', '2-somy' etc.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a
#'   file was specified.
plotProfile <- function(model, both.strands=FALSE, plot.breakpoints=FALSE, file=NULL, normalize.counts=NULL) {

    if (class(model)=='GRanges') {
        binned.data <- model
        model <- list()
        class(model) <- "aneuHMM"
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'))
        plot.profile(model, both.strands=both.strands, plot.breakpoints=FALSE, file=file)
    } else if (is(model, "GRangesList")) {
        binned.data <- model[[1]]
        model <- list()
        class(model) <- "aneuHMM"
        model$ID <- ''
        model$bins <- binned.data
        model$qualityInfo <- list(entropy=qc.entropy(binned.data$counts), spikiness=qc.spikiness(binned.data$counts), complexity=attr(binned.data, 'complexity'))
        plot.profile(model, both.strands=both.strands, plot.breakpoints=FALSE, file=file)
    } else if (class(model)=="aneuHMM") {
        plot.profile(model, both.strands=FALSE, plot.breakpoints=FALSE, file=file, normalize.counts = normalize.counts)
    } else if (class(model)=="aneuBiHMM") {
        plot.profile(model, both.strands=both.strands, plot.breakpoints=plot.breakpoints, file=file, normalize.counts = normalize.counts)
    }

}

plot.profile <- function(model, both.strands=FALSE, plot.breakpoints=TRUE, file=NULL, normalize.counts=NULL) {
    
    ## Convert to GRanges
    if (!is.null(model$bins$counts)) {
        bins <- model$bins
    } else if (!is.null(model$bincounts[[1]]$counts)) {
        bins <- model$bincounts[[1]]
    }
    ## Get breakpoint coordinates
    if (is.null(model$breakpoints) & plot.breakpoints) {
        warning("Cannot breakpoints coordinates. Please run 'getBreakpoints' first.")
        plot.breakpoints <- FALSE
    }
    if (plot.breakpoints) {
        bp.coords <- model$breakpoints
        # Set to midpoint
        start(bp.coords) <- (start(bp.coords)+end(bp.coords))/2
        end(bp.coords) <- start(bp.coords)
    }

    ## Get some variables
    num.chroms <- length(levels(seqnames(bins)))
    maxseqlength <- max(seqlengths(bins))
    tab <- table(bins$counts)
    tab <- tab[names(tab)!='0']
    if (both.strands) {
      custom.xlim <- get_rightxlim(c(bins$mcounts, bins$pcounts))
    } else {
      custom.xlim <- get_rightxlim(bins$counts)
    }

    ## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
    bins <- transCoord(bins)
    if (plot.breakpoints) {
        bp.coords <- transCoord(bp.coords)
    }

    # Plot the read counts
    dfplot <- as.data.frame(bins)
    # Set values too big for plotting to limit
    if (both.strands) {
        dfplot$mcounts <- - dfplot$mcounts    # negative minus counts
    }
    # Mean counts for CNV-state
    if (!is.null(model$segments$state)) {
        dfplot.seg <- as.data.frame(transCoord(model$segments))
        if (class(model)=="aneuHMM") {
            dfplot.seg$counts.CNV <- model$distributions[as.character(dfplot.seg$state),'mu']
        } else if (class(model)=="aneuBiHMM") {
            if (!is.null(model$distributions$both)) {
                dfplot.seg$counts.CNV <- model$distributions$both[as.character(dfplot.seg$state),'mu']
            } else {
                dfplot.seg$counts.CNV <- model$distributions$plus[as.character(dfplot.seg$state),'mu']
            }
            dfplot.seg$pcounts.CNV <- model$distributions$plus[as.character(dfplot.seg$pstate),'mu']
            dfplot.seg$mcounts.CNV <- -model$distributions$minus[as.character(dfplot.seg$mstate),'mu']
        }
    }
    # Normalize counts
    ylabstring <- 'read count'
    if (!is.null(normalize.counts)) {
        colmask <- grepl('counts', names(dfplot))
        if (class(model)=="aneuHMM") {
            nfactor <- model$distributions[normalize.counts, 'mu']
        } else if (class(model)=="aneuBiHMM") {
            nfactor <- model$distributions$plus[normalize.counts, 'mu']
        }
        dfplot[,colmask] <- dfplot[,colmask] / nfactor
        colmask <- grepl('counts', names(dfplot.seg))
        dfplot.seg[,colmask] <- dfplot.seg[,colmask] / nfactor
        custom.xlim <- custom.xlim / nfactor
        ylabstring <- 'normalized read count'
    }

    empty_theme <- theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
    ggplt <- ggplot(dfplot, aes_string(x='start.genome', y='counts'))    # data
    if (both.strands) {
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='pcounts'), position=position_jitter(width=0, height=0))    # read count
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='mcounts'), position=position_jitter(width=0, height=0))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='pcounts'))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='mcounts'))    # read count
    } else {
        ggplt <- ggplt + geom_jitter(aes_string(x='start.genome', y='counts'), position=position_jitter(width=0, height=0))    # read count
#         ggplt <- ggplt + geom_point(aes_string(x='start.genome', y='counts'))    # read count
    }
    if (!is.null(model$segments$state)) {
        if (both.strands) {
            ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='pcounts.CNV',xend='end.genome',yend='pcounts.CNV', col='pstate'), size=2)
            ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='mcounts.CNV',xend='end.genome',yend='mcounts.CNV', col='mstate'), size=2)
        } else {
            ggplt <- ggplt + geom_segment(data=dfplot.seg, mapping=aes_string(x='start.genome',y='counts.CNV',xend='end.genome',yend='counts.CNV', col='state'), size=2)
        }
    }
    # Chromosome lines
    cum.seqlengths <- cumsum(as.numeric(seqlengths(bins)))
    names(cum.seqlengths) <- seqlevels(bins)
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(bins)
    label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(bins) )
    df.chroms <- data.frame(x=c(0,cum.seqlengths))
    ggplt <- ggplt + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2)
    
    if (!is.null(model$segments$state)) {
        ggplt <- ggplt + scale_color_manual(name="state", values=stateColors(levels(dfplot.seg$state)), drop=FALSE)    # do not drop levels that are not present
    }
    if (plot.breakpoints) {
        df.bp <- as.data.frame(bp.coords)
        if (nrow(df.bp)>0) {
          statelevels <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
          suppressMessages( ggplt <- ggplt + scale_color_manual(name="state", values=c(breakpointColors(), stateColors(statelevels)), drop=FALSE) )
            ggplt <- ggplt + geom_segment(data=df.bp, aes_string(x='start.genome', xend='start.genome', color='type'), y=-1.5*custom.xlim, yend=-1.3*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'), alpha=0.5)
        }
    }
    ggplt <- ggplt + empty_theme    # no axes whatsoever
    if (both.strands) {
        ggplt <- ggplt + coord_cartesian(ylim=c(-1.5*custom.xlim,custom.xlim))    # set x- and y-limits
    } else {
        ggplt <- ggplt + coord_cartesian(ylim=c(0,custom.xlim))    # set x- and y-limits
    }
    # Get midpoints of each chromosome for xticks
    ggplt <- ggplt + scale_x_continuous(breaks=seqlengths(model$bins)/2+cum.seqlengths.0[as.character(seqlevels(model$bins))], labels=seqlevels(model$bins))
    # Quality info
    qualityInfo <- getQC(model)
    quality.string <- paste0('reads = ',round(qualityInfo$total.read.count/1e6,2),'M, complexity = ',round(qualityInfo$complexity/1e6,2),'M,  spikiness = ',round(qualityInfo$spikiness,2),',  entropy = ',round(qualityInfo$entropy,2),',  bhattacharyya = ',round(qualityInfo$bhattacharyya,2), ', num.segments = ',qualityInfo$num.segments, ', loglik = ',round(qualityInfo$loglik), ', sos = ',round(qualityInfo$sos))
    ggplt <- ggplt + ylab(ylabstring) + ggtitle(bquote(atop(.(model$ID), atop(.(quality.string),''))))
        
    if (!is.null(file)) {
        ggsave(file, ggplt, width=20, height=5)
    } else {
        return(ggplt)
    }
}
