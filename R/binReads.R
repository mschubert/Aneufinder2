binReads <- function(reads, bins) {
    UseMethod("binReads")
}

binReads.character <- function(reads, bins) {
    bin_one <- function(r, b) {
        gr = readGRanges(r)
        binReads(gr, b)
    }

    if (length(reads) == 1 && dir.exists(reads)) {
        reads = list.files(reads, "\\.(bam|bed(\\.gz)?)$", full.names=TRUE)
        message("Reading ", length(reads), " sequence files ...")
    }

    names(reads) = tools::file_path_sans_ext(basename(reads))
    binReads(as.list(reads))
}

binReads.GRanges <- function(reads, bins) {
    reads.plus <- reads[GenomicRanges::strand(reads) == '+']
    reads.minus <- reads[GenomicRanges::strand(reads) == '-']
    reads.star <- reads[GenomicRanges::strand(reads) == '*']

    scounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.star))
    mcounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.minus))
    pcounts <- suppressWarnings(GenomicRanges::countOverlaps(bins, reads.plus))

    mcols(bins) <- S4Vectors::DataFrame(counts = mcounts + pcounts + scounts,
                                        mcounts = mcounts,
                                        pcounts = pcounts)

    attr(bins, 'ID') <- attr(reads, 'ID')
    bins
}

binReads.list <- function(reads, bins) {
    lapply(reads, binReads)
}

binReads.default <- function(reads, bins) {
    stop("Do not know to to bin object of class ", sQuote(class(reads)))
}
