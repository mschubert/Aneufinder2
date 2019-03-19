#TODO: check same genome
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

binReads.list <- function(reads, bins) {
    message("Reading ", length(reads), " sequence files ...")
    lapply(reads, binReads, bins=bins)
}

binReads.default <- function(reads, bins) {
    stop("Do not know to to bin object of class ", sQuote(class(reads)))
}
