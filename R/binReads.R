binReads <- function(reads, bins) {
    reads.plus <- reads[strand(reads) == '+']
    reads.minus <- reads[strand(reads) == '-']
    reads.star <- reads[strand(reads) == '*']

    scounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.star))
    mcounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.minus))
    pcounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.plus))

    counts <- mcounts + pcounts + scounts
    countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
    colnames(countmatrix) <- c('counts','mcounts','pcounts')

    mcols(bins) <- as(countmatrix, 'DataFrame')
    attr(bins, 'ID') <- attr(reads, 'ID')
    bins
}
