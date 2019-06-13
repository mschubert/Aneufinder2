#' Genome wide heatmap of CNV-state
#'
#' Plot a genome wide heatmap of copy number variation state. This heatmap is
#' best plotted to file, because in most cases it will be too big for cleanly
#' plotting it to screen.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector
#'   with files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the
#'   same length as \code{hmms}. If \code{NULL} the IDs from the
#'   \code{\link{aneuHMM}} objects will be used.
#' @param classes A character vector with the classification of the elements on
#'   the y-axis. The vector must have the same length as \code{hmms}.
#' @param reorder.by.class If \code{TRUE}, the dendrogram will be reordered to
#'   display similar classes next to each other.
#' @param classes.color A (named) vector with colors that are used to
#'   distinguish \code{classes}. Names must correspond to the unique elements in
#'   \code{classes}.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the
#'   samples should be clustered by similarity in their CNV-state.
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @param hotspots A \code{\link{GRanges-class}} object with coordinates of
#'   genomic hotspots (see \code{\link{hotspotter}}).
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will
#'   be excluded from the computation of the clustering. This can be useful to
#'   exclude regions with artifacts.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @importFrom stats as.dendrogram
#' @importFrom ggdendro dendro_data theme_dendro
#' @importFrom S4Vectors endoapply
#' @export
#' @examples
#' ## Get results from a small-cell-lung-cancer
#' lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' lung.files <- list.files(lung.folder, full.names=TRUE)
#' ## Get results from the liver metastasis of the same patient
#' liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#' liver.files <- list.files(liver.folder, full.names=TRUE)
#' ## Plot a clustered heatmap
#' classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#' labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#' heatmapGenomewide(c(lung.files, liver.files), ylabels=labels, classes=classes,
#'                   classes.color=c('blue','red'))
heatmapGenomewide <- function(models, cluster=FALSE) {
    if ("aneuHMM" %in% class(models))
        models = list(models)

    segs <- do.call("rbind", lapply(models, function(m) {
        segs <- as.data.frame(m$segments)
        segs$ID <- m$ID
        segs
    }))

    ggplot(segs, aes(color=state)) +
        geom_segment(aes(x=start, xend=end, y=ID, yend=ID), size=8) +
        facet_grid(. ~ seqnames, scales="free_x", space="free_x") +
        scale_color_manual(values=stateColors(levels(segs$state))) +
        theme(panel.spacing = unit(0, "lines"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}
