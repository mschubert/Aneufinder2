genome <- function(x, ...) {
    UseMethod("genome")
}

genome.character <- function(x, ...) {
    if (length(x) > 1)
        stop("Can only get genome information for one object at a time")

    if (file.exists(x) && grepl("\\.bam$", x))
        genome(Rsamtools::BamFile(x), ...)
    else
        genome(GenomeInfoDb::fetchExtendedChromInfoFromUCSC(x), ...)
}

genome.GRanges <- function(x, ...) {
    seqinfo <- GenomeInfoDb::seqinfo(x)
    genome(seqinfo, ...)
}

genome.BamFile <- function(x, ...) {
    seqinfo <- GenomeInfoDb::seqinfo(x)
    seqinfo(seqinfo, ...)
}

genome.seqinfo <- function(x, chrs=NULL) {
    if (is.null(chrs)) {
        chrs = GenomeInfoDb::standardChromosomes(x)
        chrs = grep("^(chr)?[YM]$", chrs, invert=TRUE, value=TRUE)
    }
    GenomeInfoDb::keepSeqlevels(x, chrs, pruning.mode="coarse")
}

genome.default <- function(x, ...) {
    stop("Can not parse genome info from type ", sQuote(class(x)))
}
