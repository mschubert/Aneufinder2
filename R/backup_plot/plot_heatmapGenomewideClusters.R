
#' Plot heatmaps for quality control
#' 
#' This function is a convenient wrapper to call
#' \code{\link{heatmapGenomewide}} for all clusters after calling
#' \code{\link{clusterByQuality}} and plot the heatmaps into one pdf for
#' efficient comparison.
#' 
#' @param cl The return value of \code{\link{clusterByQuality}}.
#' @param cutree The return value of \code{\link[stats]{cutree}}, where the
#'   names correspond to the filenames to be loaded.
#' @param file A character specifying the output file.
#' @param ... Further parameters passed on to \code{\link{heatmapGenomewide}}.
#' @return A \code{\link[cowplot]{cowplot}} object or \code{NULL} if a file was specified.
#' @export
#' @examples
#' ## Get a list of HMMs and cluster them
#' folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' files <- list.files(folder, full.names=TRUE)
#' cl <- clusterByQuality(files, G=5)
#' heatmapGenomewideClusters(cl=cl)
#' 
#' ## Plot sub-clones of the largest cluster
#' largest.cluster <- which.max(sapply(cl$classification, length))
#' files <- cl$classification[[largest.cluster]]
#' clust <- clusterHMMs(files)
#' groups <- cutree(tree = clust$hclust, k = 5)
#' heatmapGenomewideClusters(cutree = groups, cluster = FALSE)
heatmapGenomewideClusters <- function(cl=NULL, cutree=NULL, file=NULL, ...) {
  
    ## Check user input ##
    if (is.null(cl) & is.null(cutree)) {
        stop("Please specify either 'cl' or 'cutree'.")
    }
    if (!is.null(cl) & !is.null(cutree)) {
        stop("Please specify either 'cl' or 'cutree', not both.")
    }
  
    if (!is.null(cl)) {
        filelist <- cl$classification
    } else if (!is.null(cutree)) {
        filelist <- split(names(cutree), cutree)
    }
    ## Get the plot dimensions ##
    ptm <- startTimedMessage("Calculating plot dimensions ...")
    hmm <- loadFromFiles(filelist[[1]][1])[[1]]
      width.heatmap <- sum(as.numeric(seqlengths(hmm$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
      height <- max(length(unlist(filelist)) * 0.5, 2)
      width.dendro <- 20
      width <- width.heatmap + width.dendro
    stopTimedMessage(ptm)
  
    ## Make the plots ##
    ggplts <- list()
    for (i1 in 1:length(filelist)) {
        message("Cluster ", i1, " ...")
        ggplts[[i1]] <- heatmapGenomewide(filelist[[i1]], ...)
    }
    cowplt <- cowplot::plot_grid(plotlist = ggplts, align='v', ncol=1, rel_heights = sapply(filelist, function(x) { max(length(x), 4) }))
    if (is.null(file)) {
        return(cowplt)
    } else {
        ggsave(cowplt, filename = file, width=width, height=height, units='cm', limitsize = FALSE)
    }
    
}
