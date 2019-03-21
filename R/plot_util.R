#' @import ggplot2
#' @importFrom cowplot plot_grid draw_label
#' @import reshape2
NULL

#' \pkg{AneuFinder} color scheme
#'
#' Get the color schemes that are used in the AneuFinder plots.
#'
#' @return A character vector with colors.
#' @name colors
#' @aliases stateColors strandColors
NULL

#' Plotting function for binned read counts
#'
#' Make plots for binned read counts from \code{\link{binned.data}}.
#'
#' @param x A \code{\link{GRanges-class}} object with binned read counts.
#' @param type Type of the plot, one of \code{c('profile', 'histogram',
#'   'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with read counts.}
#'   \item{\code{histogram}}{A histogram of read counts.}
#'   \item{\code{profile}}{An profile with read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRanges
#' @export
plot.GRanges <- function(x, type='profile', ...) {
    if (type == 'karyogram' | type==3) {
        plotKaryogram(x, ...)
    } else if (type == 'histogram' | type==2) {
        plotHistogram(x, ...)
    } else if (type == 'profile' | type==1) {
        plotProfile(x, ...)
    }
}

#' Plotting function for binned read counts (list)
#'
#' Make plots for binned read counts (list) from \code{\link{binned.data}}.
#'
#' @param x A \code{\link{GRangesList}} object with binned read counts.
#' @param type Type of the plot, one of \code{c('profile', 'histogram',
#'   'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with read counts.}
#'   \item{\code{histogram}}{A histogram of read counts.}
#'   \item{\code{profile}}{An profile with read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRangesList
#' @export
plot.GRangesList <- function(x, type='profile', ...) {
    if (type == 'karyogram' | type==3) {
        plotKaryogram(x, ...)
    } else if (type == 'histogram' | type==2) {
        plotHistogram(x, ...)
    } else if (type == 'profile' | type==1) {
        plotProfile(x, ...)
    }
}

#' State colors
#'
#' @describeIn colors Colors that are used for the states.
#' @param states A character vector with states whose color should be returned.
#' @export
#' @examples
#' ## Make a nice pie chart with the AneuFinder state color scheme
#' statecolors <- stateColors()
#' pie(rep(1,length(statecolors)), labels=names(statecolors), col=statecolors)
stateColors <- function(states=c('zero-inflation', paste0(0:10, '-somy'), 'total')) {
    state.colors <- c("zero-inflation"="gray90", "0-somy"="gray90","1-somy"="darkorchid3","2-somy"="springgreen2","3-somy"="red3","4-somy"="gold2","5-somy"="navy","6-somy"="lemonchiffon","7-somy"="dodgerblue","8-somy"="chartreuse4","9-somy"="lightcoral","10-somy"="aquamarine2","total"="black")
    states.with.color <- intersect(states, names(state.colors))
    cols <- rep('black', length(states))
    names(cols) <- states
    cols[states.with.color] <- state.colors[states.with.color]
    return(cols)
}

#' Strand colors
#'
#' @describeIn colors Colors that are used to distinguish strands.
#' @param strands A character vector with strands whose color should be returned. Any combination of \code{c('+','-','*')}.
#' @export
#' @examples
#' ## Make a nice pie chart with the AneuFinder strand color scheme
#' strandcolors <- strandColors()
#' pie(rep(1,length(strandcolors)), labels=names(strandcolors), col=strandcolors)
strandColors <- function(strands=c('+','-')) {
    strand.colors <- c('+'="#678B8B", '-'="#F3A561", '*'="#000000")
    strands.with.color <- intersect(strands, names(strand.colors))
    cols <- rep('black', length(strands))
    names(cols) <- strands
    cols[strands.with.color] <- strand.colors[strands.with.color]
    return(cols)
}

#' Breakpoint colors
#'
#' @describeIn colors Colors that are used for breakpoint types.
#' @param breaktypes A character vector with breakpoint types whose color
#' should be returned. Any combination of
#' \code{c('CNB','SCE','CNB+SCE','other')}.
#' @export
#' @examples
#' ## Make a nice pie chart with the AneuFinder breakpoint-type color scheme
#' breakpointcolors <- breakpointColors()
#' pie(rep(1,length(breakpointcolors)), labels=names(breakpointcolors), col=breakpointcolors)
breakpointColors <- function(breaktypes=c('CNB','SCE','CNB+SCE','other')) {
    break.colors <- c('CNB'="dodgerblue4", 'SCE'="tomato3", 'CNB+SCE'="orchid4", 'other'='gray30')
    breaks.with.color <- intersect(breaktypes, names(break.colors))
    cols <- rep('gray30', length(breaktypes))
    names(cols) <- breaktypes
    cols[breaks.with.color] <- break.colors[breaks.with.color]
    return(cols)
}

#' Plotting function for saved \pkg{\link{AneuFinder}} objects
#'
#' Convenience function that loads and plots a \pkg{\link{AneuFinder}} object in one step.
#'
#' @param x A filename that contains either \code{\link{binned.data}} or a
#'   \code{\link{aneuHMM}}.
#' @param ... Additional arguments.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot character
#' @importFrom graphics plot
#' @export
plot.character <- function(x, ...) {
    x <- get(load(x))
    graphics::plot(x, ...)
}

get_rightxlim <- function(counts) {
#     rightxlim1 <- median(counts[counts>0])*7
#     tab <- table(counts)
#     tab <- tab[names(tab)!='0']
#     breaks <- as.numeric(names(tab))
#     rightxlim2 <- breaks[tab<=5 & breaks>median(counts)*2][1]
#     rightxlim <- min(rightxlim1,rightxlim2, na.rm=TRUE)
    rightxlim <- stats::quantile(counts, 0.999)
    if (length(rightxlim)==0 | is.na(rightxlim) | is.infinite(rightxlim)) {
        rightxlim <- 1
    }
    return(rightxlim)
}

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the
#' \code{\link{GRanges-class}} object. This is useful for making genomewide
#' plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} with two additional metadata
#' columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
    cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(gr)
    gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    return(gr)
}

#' Plotting function for \code{\link{aneuHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuHMM}} objects.
#'
#' @param x An \code{\link{aneuHMM}} object.
#' @param type Type of the plot, one of \code{c('profile', 'histogram',
#'   'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{An profile with read counts and CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuHMM
#' @export
plot.aneuHMM <- function(x, type='profile', ...) {
    if (type == 'karyogram' | type==3) {
        plotKaryogram(x, ...)
    } else if (type == 'histogram' | type==2) {
        plotHistogram(x, ...)
    } else if (type == 'profile' | type==1) {
        plotProfile(x, ...)
    }
}

#' Plotting function for \code{\link{aneuBiHMM}} objects
#'
#' Make different types of plots for \code{\link{aneuBiHMM}} objects.
#'
#' @param x An \code{\link{aneuBiHMM}} object.
#' @param type Type of the plot, one of \code{c('profile', 'histogram',
#'   'karyogram')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{profile}}{An profile with read counts and CNV-state.}
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{A karyogram-like chromosome overview with CNV-state.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot aneuBiHMM
#' @export
plot.aneuBiHMM <- function(x, type='profile', ...) {
    if (type == 'karyogram' | type==3) {
        args <- names(list(...))
        if ('both.strands' %in% args) {
            plotKaryogram(x, ...)
        } else {
            plotKaryogram(x, both.strands=TRUE, ...)
        }
    } else if (type == 'histogram' | type==2) {
        plotBivariateHistograms(x, ...)
    } else if (type == 'profile' | type==1) {
        args <- names(list(...))
        if ('both.strands' %in% args) {
            plotProfile(x, ...)
        } else {
            plotProfile(x, both.strands=TRUE, ...)
        }
    }
}
