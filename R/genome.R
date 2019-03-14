genome <- function(x, ...) {
    UseMethod("genome")
}

genome.character <- function(x, ...) {
    if (length(x) > 1)
        stop("Can only get genome information for one object at a time")
    if (file.exists(x) && grepl("\\.bam$", x))
        genome(Rsamtools::BamFile(x), ...)

    #TODO: get seqinfo object from assembly identifier
}

genome.GRanges <- function(x, ...) {
    seqinfo <- GenomeInfoDB::seqinfo(x)
    genome(seqinfo, ...)
}

genome.BamFile <- function(x, ...) {
    seqinfo <- GenomeInfoDB::seqinfo(x)
    seqinfo(seqinfo, ...)
}

genome.seqinfo <- function(x, ...) {
    #TODO: remove nonstandard chromosomes by default
}

genome.default <- function(x, ...) {
    stop("Can not parse genome info from type ", sQuote(class(x)))
}
