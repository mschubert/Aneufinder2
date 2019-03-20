#' Convert aligned reads from various file formats into read counts in
#' equidistant bins
#'
#' Convert aligned reads in .bam or .bed(.gz) format into read counts in
#' equidistant windows.
#'
#' Convert aligned reads from .bam or .bed(.gz) files into read counts in
#' equidistant windows (bins). This function uses
#' \code{GenomicRanges::countOverlaps} to calculate the read counts.
#'
#' @param file A file with aligned reads. Alternatively a
#'   \code{\link{GRanges-class}} with aligned reads.
#' @param ID An identifier that will be used to identify the file throughout
#'   the workflow and in plotting.
#' @inheritParams readGRanges
#' @param outputfolder.binned Folder to which the binned data will be saved. If
#'   the specified folder does not exist, it will be created.
#' @param binsizes An integer vector with bin sizes. If more than one value is
#'   given, output files will be produced for each bin size.
#' @param stepsizes A vector of step sizes the same length as \code{binsizes}.
#'   Only used for \code{method="HMM"}.
#' @param bins A named \code{list} with \code{\link{GRanges-class}} containing
#'   precalculated bins produced by \code{\link{fixedWidthBins}} or
#'   \code{\link{variableWidthBins}}. Names must correspond to the binsize.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin
#'   size will be selected accordingly. Output files are produced for each value.
#' @param reads.per.step Approximate number of desired reads per step.
#' @param variable.width.reference A BAM file that is used as reference to
#'   produce variable width bins. See \code{\link{variableWidthBins}} for
#'   details.
#' @param chromosomes If only a subset of the chromosomes should be binned,
#'   specify them here.
#' @param save.as.RData If set to \code{FALSE}, no output file will be written.
#'   Instead, a \link{GenomicRanges} object containing the binned data will be
#'   returned. Only the first binsize will be processed in this case.
#' @param calc.complexity A logical indicating whether or not to estimate
#'   library complexity.
#' @param call The \code{match.call()} of the parent function.
#' @param reads.store If \code{TRUE} processed read fragments will be saved to
#'   file. Reads are processed according to \code{min.mapq} and
#'   \code{remove.duplicate.reads}. Paired end reads are coerced to single end
#'   fragments. Will be ignored if \code{use.bamsignals=TRUE}.
#' @param outputfolder.reads Folder to which the read fragments will be saved.
#'   If the specified folder does not exist, it will be created.
#' @param reads.return If \code{TRUE} no binning is done and instead, read
#'   fragments from the input file are returned in \code{\link{GRanges-class}}
#'   format.
#' @param reads.overwrite Whether or not an existing file with read fragments
#'   should be overwritten.
#' @param reads.only If \code{TRUE} only read fragments are stored and/or
#'   returned and no binning is done.
#' @param use.bamsignals If \code{TRUE} the \pkg{\link[bamsignals]{bamsignals}}
#'   package will be used for binning. This gives a tremendous performance
#'   increase for the binning step. \code{reads.store} and \code{calc.complexity}
#'   will be set to \code{FALSE} in this case.
#' @return The function produces a \code{list()} of \code{\link{GRanges-class}}
#'   or \code{\link{GRangesList}} objects with meta data columns 'counts',
#'   'mcounts', 'pcounts' that contain the total, minus and plus read count. This
#'   binned data will be either written to file (\code{save.as.RData=FALSE}) or
#'   given as return value (\code{save.as.RData=FALSE}).
#' @seealso binning
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom bamsignals bamCount
#' @export
#'
#' @examples
#' ## Get an example BED file with single-cell-sequencing reads
#' bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#' ## Bin the BED file into bin size 1Mb
#' binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                    chromosomes=c(1:19,'X','Y'))
#' print(binned)
binReads <- function(reads, bins) {
    UseMethod("binReads")
}

binReads.character <- function(reads, bins) {
    if (length(reads) == 1 && dir.exists(reads))
        reads = list.files(reads, "\\.(bam|bed(\\.gz)?)$", full.names=TRUE)
    if (length(reads) == 0)
        stop("No BAM or BED files supplied or in directory")
    names(reads) = tools::file_path_sans_ext(basename(reads))

    if (length(reads) > 1)
        binReads(as.list(reads), bins)
    else
        binReads(readGRanges(reads), bins)
}

binReads.GRanges <- function(reads, bins) {
    strands <- split(reads, GenomicRanges::strand(reads))
    clist <- lapply(strands, function(r)
                    suppressWarnings(GenomicRanges::countOverlaps(bins, r)))

    bins$counts <- Reduce(`+`, clist)
    bins$mcounts <- clist$`-`
    bins$pcounts <- clist$`+`

    attr(bins, 'ID') <- attr(reads, 'ID')
    bins
}

#TODO: check same genome
binReads.list <- function(reads, bins) {
    message("Reading ", length(reads), " sequence files ...")
    lapply(reads, binReads, bins=bins)
}

binReads.default <- function(reads, bins) {
    stop("Do not know to to bin object of class ", sQuote(class(reads)))
}
