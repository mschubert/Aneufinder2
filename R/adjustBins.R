#' Make variable-width bins
#' 
#' Make variable-width bins based on a reference BAM file. This can be a
#' simulated file (produced by \code{\link{simulateReads}} and aligned with
#' your favourite aligner) or a real reference.
#' 
#' Variable-width bins are produced by first binning the reference BAM file
#' with fixed-width bins and selecting the desired number of reads per bin as
#' the (non-zero) maximum of the histogram. A new set of bins is then generated
#' such that every bin contains the desired number of reads.
#' 
#' @param reads A \code{\link{GRanges-class}} with reads. See
#'   \code{\link{bam2GRanges}} and \code{\link{bed2GRanges}}.
#' @param binsizes A vector with binsizes. Resulting bins will be close to the
#'   specified binsizes.
#' @param stepsizes A vector of step sizes in base pairs, the same length as
#'   \code{binsizes}.
#' @param chromosomes A subset of chromosomes for which the bins are generated.
#' @return A \code{list()} of \code{\link{GRanges-class}} objects with
#'   variable-width bins. If \code{stepsizes} is specified, a \code{list()} of
#'   \code{\link{GRangesList}} objects with one entry per step.
#' @author Aaron Taudt
#' @importFrom S4Vectors endoapply
#' @export
#'
#' @examples
#' ## Get an example BED file with single-cell-sequencing reads
#' bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#' ## Read the file into a GRanges object
#' reads <- bed2GRanges(bedfile, assembly='mm10', chromosomes=c(1:19,'X','Y'),
#'                      min.mapq=10, remove.duplicate.reads=TRUE)
#' ## Make variable-width bins of size 500kb and 1Mb
#' bins <- variableWidthBins(reads, binsizes=c(5e5,1e6))
#' ## Plot the distribution of binsizes
#' hist(width(bins[['binsize_1e+06']]), breaks=50)
adjustBins <- function(bins, reads) {
    bins <- binReads(reads, bins)
    strand(reads) <- '*'
    reads <- sort(reads)

    binsize <- attr(bins, 'binsize')
    stepsize <- attr(bins, 'stepsize')
    numsteps <- binsize / stepsize

    mediancount <- as.integer(median(bins$counts[bins$counts>0]))
    mc.perstep <- unique(as.integer((1:numsteps) * mediancount / numsteps))

    # Determine consecutive bins for each shift.
    for (istep in 1:length(mc.perstep)) {
        subreads <- GRangesList()
        skipped.chroms <- character()

        # Pick only every mediancount read for each chromosome separately.
        for (chrom in unique(seqnames(bins))) {
            reads.chr <- reads[seqnames(reads)==chrom]
            if (length(reads.chr) >= mediancount){
                idx <- seq(mc.perstep[istep], length(reads.chr), by=mediancount)
                subreads[[chrom]] <- reads.chr[idx]
            } else {
                skipped.chroms[chrom] <- chrom
            }
        }

        subreads <- unlist(subreads, use.names=FALSE)
        # Adjust length of reads to get consecutive bins.
        subreads <- resize(subreads, width=1)
        # Make new bins. Gaps until seqlengths-1 because we have to add 1 later
        # to get consecutive bins
        bins <- gaps(subreads, start=1L, end=seqlengths(subreads)-1L)
        bins <- bins[strand(bins)=='*']
        end(bins) <- end(bins) + 1
        bins.split <- split(bins, seqnames(bins))

        # Remove first bin for smaller than median count steps (overlapping bins).
        if (mc.perstep[istep] < mediancount) {
            bins.split <- endoapply(bins.split, function(x) { x[-1] })
        }

        # We don't want incomplete bins (number of reads < median count) at the
        # end of each chromosome.
        bins.split <- endoapply(bins.split, function(x) { x[-length(x)] })
        bins <- unlist(bins.split, use.names=FALSE)
        bins <- bins[!seqnames(bins) %in% skipped.chroms]
        bins <- keepSeqlevels(bins, setdiff(seqlevels(bins), skipped.chroms))
    }

    bins
}
