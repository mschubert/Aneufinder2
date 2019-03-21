#' Perform a PCA for copy number profiles
#'
#' Perform a PCA for copy number profiles in \code{\link{aneuHMM}} objects.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector
#'   with files that contain such objects.
#' @param PC1 Integer specifying the first of the principal components to plot.
#' @param PC2 Integer specifying the second of the principal components to plot.
#' @param colorBy A character vector of the same length as \code{hmms} which is
#'   used to color the points in the plot.
#' @param plot Set to \code{FALSE} if you want to return the data.frame that is
#'   used for plotting instead of the plot.
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will
#'   be excluded from the computation of the PCA. This can be useful to exclude
#'   regions with artifacts.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or a data.frame if
#'   \code{plot=FALSE}.
#' @export
#' @examples
#' ## Get results from a small-cell-lung-cancer
#' lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' lung.files <- list.files(lung.folder, full.names=TRUE)
#' ## Get results from the liver metastasis of the same patient
#' liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#' liver.files <- list.files(liver.folder, full.names=TRUE)
#' ## Plot the PCA
#' classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#' labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#' plot_pca(c(lung.files, liver.files), colorBy=classes, PC1=2, PC2=4)
plot_pca <- function(hmms, PC1=1, PC2=2, colorBy=NULL, plot=TRUE,
                     exclude.regions=NULL) {
  
    hmms <- loadFromFiles(hmms, check.class = c("aneuHMM", "aneuBiHMM"))
    copy.numbers <- sapply(hmms, function(hmm) { hmm$bins$copy.number })
    
    ### Exclude regions ###
    if (!is.null(exclude.regions)) {
        ind <- findOverlaps(hmms[[1]]$bins, exclude.regions)@from
            copy.numbers <- copy.numbers[-ind,]
    }
    
    ## PCA
    Y <- apply(log(copy.numbers+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
    s <- svd(Y)
    percent <- s$d^2/sum(s$d^2)*100
    labs <- sapply(seq_along(percent), function(i) {
        paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
    })
    df <- data.frame(PC1 = s$u[,PC1], PC2 = s$u[,PC2])
    if (!plot) {
        colnames(df) <- labs[c(PC1, PC2)]
        rownames(df) <- colnames(copy.numbers)
        return(df)
    } else {
        if (!is.null(colorBy)) {
            df$color <- colorBy
            ggplt <- ggplot(df) + geom_point(aes_string(x='PC1', y='PC2', color='color'))
            ggplt <- ggplt + scale_color_manual(values=getDistinctColors(unique(colorBy)))
        } else {
            ggplt <- ggplot(df) + geom_point(aes_string(x='PC1', y='PC2'))
        }
        ggplt <- ggplt + xlab(labs[PC1]) + ylab(labs[PC2])
        return(ggplt)
    }
}
