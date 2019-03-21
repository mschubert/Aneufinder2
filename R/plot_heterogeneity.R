#' Heterogeneity vs. Aneuploidy
#' 
#' Make heterogeneity vs. aneuploidy plots using individual chromosomes as
#' datapoints.
#' 
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector
#'   with files that contain such objects.
#' @param hmms.list Alternative input for a faceted plot. A named list() of
#'   lists of \code{\link{aneuHMM}} objects. Alternatively a named list() of
#'   character vectors with files that contain \code{\link{aneuHMM}} objects.
#'   List names serve as facets for plotting. If specified,
#'   \code{normalChromosomeNumbers} is assumed to be a list() of vectors (or
#'   matrices) instead of a vector (or matrix).
#' @param normalChromosomeNumbers A named integer vector or matrix with
#'   physiological copy numbers, where each element (vector) or column (matrix)
#'   corresponds to a chromosome. This is useful to specify male or female
#'   samples, e.g. \code{c('X'=2)} for female samples or \code{c('X'=1,'Y'=1)}
#'   for male samples. Specify a vector if all your \code{hmms} have the same
#'   physiological copy numbers. Specify a matrix if your \code{hmms} have
#'   different physiological copy numbers (e.g. a mix of male and female
#'   samples). If not specified otherwise, '2' will be assumed for all
#'   chromosomes. If you have specified \code{hmms.list} instead of \code{hmms},
#'   \code{normalChromosomeNumbers} is assumed to be a list() of vectors (or
#'   matrices), with one vector (or matrix) for each element in \code{hmms.list}.
#' @param plot A logical indicating whether to plot or to return the underlying data.frame.
#' @inheritParams karyotypeMeasures
#' @return A \code{\link[ggplot2]{ggplot}} object or a data.frame if \code{plot=FALSE}.
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples 
#' ### Example 1: A faceted plot of lung and liver cells ###
#' ## Get results from a small-cell-lung-cancer
#' lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' lung.files <- list.files(lung.folder, full.names=TRUE)
#' ## Get results from the liver metastasis of the same patient
#' liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#' liver.files <- list.files(liver.folder, full.names=TRUE)
#' ## Make heterogeneity plots
#' plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files))
#'
#' ### Example 2: Plot a mixture sample of male and female cells ###
#' ## Get results from a small-cell-lung-cancer
#' folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' files <- list.files(lung.folder, full.names=TRUE)
#' ## Construct a matrix with physiological copy numbers for a mix of 48 male and 48 female samples
#' normal.chrom.numbers <- matrix(2, nrow=96, ncol=24,
#'                                dimnames=list(sample=c(paste('male', 1:48), paste('female', 49:96)),
#'                                              chromosome=c(1:22,'X','Y')))
#' normal.chrom.numbers[1:48,c('X','Y')] <- 1
#' normal.chrom.numbers[49:96,c('Y')] <- 0
#' head(normal.chrom.numbers)
#' ## Make heterogeneity plots
#' plotHeterogeneity(hmms = files, normalChromosomeNumbers = normal.chrom.numbers)
#'
#' ### Example 3: A faceted plot of male lung and female liver cells ###
#' ## Get results from a small-cell-lung-cancer
#' lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' lung.files <- list.files(lung.folder, full.names=TRUE)
#' ## Specify the physiological copy numbers
#' chrom.numbers.lung <- c(rep(2, 22), 1, 1)
#' names(chrom.numbers.lung) <- c(1:22, 'X', 'Y')
#' print(chrom.numbers.lung)
#' ## Get results from the liver metastasis of the same patient
#' liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#' liver.files <- list.files(liver.folder, full.names=TRUE)
#' ## Specify the physiological copy numbers
#' chrom.numbers.liver <- c(rep(2, 22), 2, 0)
#' names(chrom.numbers.liver) <- c(1:22, 'X', 'Y')
#' print(chrom.numbers.liver)
#' ## Make heterogeneity plots
#' plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files),
#'                   normalChromosomeNumbers = list(chrom.numbers.lung, chrom.numbers.liver))
#'
#' ### Example 4 ###
#' ## Exclude artifact regions with high variance
#' consensus <- consensusSegments(c(lung.files, liver.files))
#' variance <- apply(consensus$copy.number, 1, var)
#' exclude.regions <- consensus[variance > quantile(variance, 0.999)]
#' ## Make heterogeneity plots
#' plotHeterogeneity(hmms.list = list(lung=lung.files, liver=liver.files),
#'                   exclude.regions=exclude.regions)
plotHeterogeneity <- function(hmms, hmms.list=NULL,
                              normalChromosomeNumbers=NULL, plot=TRUE,
                              regions=NULL, exclude.regions=NULL) {

    if (!is.null(hmms.list)) {
        if (!is.null(normalChromosomeNumbers)) {
            if (!class(normalChromosomeNumbers) == 'list') {
                stop("Argument 'normalChromosomeNumbers' has to be a list with one entry (vector or matrix) for each entry in 'hmms.list'.")
            }
        }
    }
    if (is.null(hmms.list)) {
        hmms <- loadFromFiles(hmms, check.class="aneuHMM")
        ## Karyotype measures
        kmeasures <- karyotypeMeasures(hmms, normalChromosomeNumbers = normalChromosomeNumbers, regions = regions, exclude.regions = exclude.regions)
        rownames(kmeasures$genomewide) <- 'all'
        kmeasures <- rbind(kmeasures$genomewide, kmeasures$per.chromosome)
        kmeasures$chromosome <- rownames(kmeasures)
        rownames(kmeasures) <- NULL
        
        if (plot) {
            ## Plot with ggrepel
            ggplt <- ggplot(data=kmeasures, mapping=aes_string(x='Aneuploidy', y='Heterogeneity')) + geom_point()
            ggplt <- ggplt + geom_text_repel(aes_string(label='chromosome'))
            return(ggplt)
        } else {
            return(kmeasures)
        }
    } else {
        kmeasures.all <- list()
        for (i1 in 1:length(hmms.list)) {
            hmms <- hmms.list[[i1]]
            samplename <- names(hmms.list)[i1]
            hmms <- loadFromFiles(hmms, check.class="aneuHMM")
            ## Karyotype measures
            kmeasures <- karyotypeMeasures(hmms, normalChromosomeNumbers = normalChromosomeNumbers[[i1]], regions = regions, exclude.regions = exclude.regions)
            rownames(kmeasures$genomewide) <- 'all'
            kmeasures <- rbind(kmeasures$genomewide, kmeasures$per.chromosome)
            kmeasures$chromosome <- rownames(kmeasures)
            kmeasures$sample <- samplename
            kmeasures.all[[i1]] <- kmeasures
        }
        kmeasures.all <- do.call(rbind, kmeasures.all)
        rownames(kmeasures.all) <- NULL
        
        if (plot) {
            ## Plot with ggrepel
            ggplt <- ggplot(data=kmeasures.all, mapping=aes_string(x='Aneuploidy', y='Heterogeneity')) + geom_point()
            ggplt <- ggplt + geom_text_repel(aes_string(label='chromosome'))
            ggplt <- ggplt + facet_wrap(~sample)
            return(ggplt)
        } else {
            return(kmeasures.all)
        }
    }
    
}
