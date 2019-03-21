#' Plot aneuploidy state
#'
#' Plot a heatmap of aneuploidy state for multiple samples. Samples can be
#' clustered and the output can be returned as data.frame.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector
#'   with files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the
#'   same length as \code{hmms}. If \code{NULL} the IDs from the
#'   \code{\link{aneuHMM}} objects will be used.
#' @param cluster If \code{TRUE}, the samples will be clustered by similarity
#'   in their CNV-state.
#' @param as.data.frame If \code{TRUE}, instead of a plot, a data.frame with
#'   the aneuploidy state for each sample will be returned.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or a data.frame,
#'   depending on option \code{as.data.frame}.
#' @author Aaron Taudt
#' @importFrom stats aggregate dist hclust
#' @export
#' @examples
#' ## Get results from a small-cell-lung-cancer
#' folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' files <- list.files(folder, full.names=TRUE)
#' ## Plot the ploidy state per chromosome
#' heatmapAneuploidies(files, cluster=FALSE)
#' ## Return the ploidy state as data.frame
#' df <- heatmapAneuploidies(files, cluster=FALSE, as.data.frame=TRUE)
#' head(df)
heatmapAneuploidies <- function(hmms, ylabels=NULL, cluster=TRUE, as.data.frame=FALSE) {

    ## Check user input
    if (!is.null(ylabels)) {
        if (length(ylabels) != length(hmms)) {
            stop("length(ylabels) must equal length(hmms)")
        }
    }
  if (length(hmms) == 1 & cluster==TRUE) {
    cluster <- FALSE
    warning("Cannot do clustering because only one object was given.")
  }

    ## Load the files
    hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
    levels.state <- unique(unlist(lapply(hmms, function(hmm) { levels(hmm$bins$state) })))
    
    ## Assign new IDs
    if (!is.null(ylabels)) {
        for (i1 in 1:length(hmms)) {
            hmms[[i1]]$ID <- ylabels[i1]
        }
    }

    ## Transform to GRanges in reduced representation
    grlred <- GRangesList()
    for (hmm in hmms) {
        if (!is.null(hmm$segments)) {
            grlred[[hmm$ID]] <- hmm$segments
        }
    }
    
    ## Find the most frequent state (mfs) for each chromosome and sample
    ptm <- startTimedMessage("finding most frequent state for each sample and chromosome ...")
    grl.per.chrom <- lapply(grlred, function(x) { split(x, seqnames(x)) })
    mfs.samples <- list()
    for (i1 in 1:length(grlred)) {
        mfs.samples[[names(grlred)[i1]]] <- list()
        for (i2 in 1:length(grl.per.chrom[[i1]])) {
          x <- grl.per.chrom[[i1]][[i2]]
      if (length(x)>0) {
        tab <- stats::aggregate(width(x), by=list(state=x$state), FUN="sum")
            mfs.samples[[names(grlred)[i1]]][[i2]] <- tab$state[which.max(tab$x)]
      } else {
            mfs.samples[[names(grlred)[i1]]][[i2]] <- "0-somy"
      }
        }
        attr(mfs.samples[[names(grlred)[i1]]], "varname") <- 'chromosome'
    }
    attr(mfs.samples, "varname") <- 'sample'
    stopTimedMessage(ptm)

    ## Transform to data.frame
    # Long format
    df <- reshape2::melt(mfs.samples, value.name='state')
    df$state <- factor(df$state, levels=levels.state)
    df$sample <- factor(df$sample, levels=unique(df$sample))
    df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
    # Wide format
    df.wide <- reshape2::dcast(df, sample ~ chromosome, value.var='state', factorsAsStrings=FALSE)
    # Correct strings to factors
    for (col in 2:ncol(df.wide)) {
        df.wide[,col] <- factor(df.wide[,col], levels=levels.state)
    }

    ## Cluster the samples by chromosome state
    if (cluster) {
        # Cluster
        hc <- stats::hclust(stats::dist(data.matrix(df.wide[-1])))
        # Reorder samples in mfs list
        mfs.samples.clustered <- mfs.samples[hc$order]
        attr(mfs.samples.clustered, "varname") <- 'sample'
        df <- reshape2::melt(mfs.samples.clustered, value.name='state')
        df$state <- factor(df$state, levels=levels.state)
        df$sample <- factor(df$sample, levels=unique(df$sample))
        df$chromosome <- factor(df$chromosome, levels=unique(df$chromosome))
    }

    ## Plot to heatmap
    if (as.data.frame) {
        df.table <- df.wide
        for (i1 in 2:ncol(df.table)) {
            df.table[,i1] <- initializeStates(levels.state)$multiplicity[df.wide[,i1]]
        }
        return(df.table)
    } else {
        ## Reorder state levels for the legend
        df$state <- factor(df$state, levels=names(sort(initializeStates(levels(df$state))$multiplicity)))
        ggplt <- ggplot(df) + geom_tile(aes_string(x='chromosome', y='sample', fill='state'), col='black') + theme_bw() + scale_fill_manual(values=stateColors(levels(df$state)))
        return(ggplt)
    }
}
