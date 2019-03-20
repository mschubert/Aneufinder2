#' GC correction
#'
#' Correct a list of \code{\link{binned.data}} by GC content.
#'
#' Two methods are available for GC correction:  Option
#' \code{method='quadratic'} uses the method described in the Supplementary of
#' \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit
#' to adjust the read count.
#' 
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects
#'   or a list of filenames containing such objects.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence
#'   that is used for the GC correction.
#' @param same.binsize If \code{TRUE} the GC content will only be calculated
#'   once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects
#'   describe the same genome at the same binsize and stepsize.
#' @param method One of \code{c('quadratic', 'loess')}. Option
#'   \code{method='quadratic'} uses the method described in the Supplementary of
#'   \code{citation("AneuFinder")}. Option \code{method='loess'} uses a loess fit
#'   to adjust the read count.
#' @param return.plot Set to \code{TRUE} if plots should be returned for visual
#'   assessment of the GC correction.
#' @param bins A \code{\link{binned.data}} object with meta-data column 'GC'.
#'   If this is specified, \code{GC.BSgenome} is ignored. Beware, no format
#'   checking is done.
#' @return A \code{list()} with \code{\link{binned.data}} objects with adjusted
#'   read counts. Alternatively a \code{list()} with
#'   \code{\link[ggplot2]{ggplot}} objects if \code{return.plot=TRUE}.
#' @author Aaron Taudt
#' @importFrom Biostrings Views alphabetFrequency
#' @importFrom stats lm predict loess
#' @importFrom reshape2 melt
#' @export
#' @examples
#' ## Get a BED file, bin it and run GC correction
#' bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#' binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                    chromosomes=c(1:19,'X','Y'))
#' plot(binned[[1]], type=1)
#' if (require(BSgenome.Mmusculus.UCSC.mm10)) {
#'   binned.GC <- correctGC(list(binned[[1]]), GC.BSgenome=BSgenome.Mmusculus.UCSC.mm10)
#'   plot(binned.GC[[1]], type=1)
#' }
correctGC <- function(reads, method='loess') {
    UseMethod("correctGC")
}

# check for same genome if passed a list?
# should this addGCcontent if not already present?
correctGC.GRanges <- function(reads, method='loess') {
    if (method == 'loess') {
        # check: use predict() here instead?
        fit <- stats::loess(counts ~ GC.content, data=as.data.frame(reads))
        correction <- median(reads$counts) / fit$fitted
        reads$counts <- as.integer(round(reads$counts * correction))
        reads$mcounts <- as.integer(round(reads$mcounts * correction))
        reads$pcounts <- as.integer(round(reads$pcounts * correction))
    }

    reads
}

correctGC.list <- function(reads, method='loess') {
    lapply(reads, correctGC, method=method)
}

correctGC.default <- function(reads, method='loess') {
    stop("Do not know how to correct GC for type ", sQuote(class(reads)))
}
