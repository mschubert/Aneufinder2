#TODO: remove nonstandard chromosomes by default

genome <- function(x, ...) {
    UseMethod("genome")
}

genome.character <- function(x, ...) {
    #TODO: get seqinfo object from assembly identifier
}

genome.GRanges <- function(x, ...) {
    GenomeInfoDB::seqinfo(x, ...)
}

genome.default <- function(x, ...) {
    stop("Can not parse genome info from type ", sQuote(class(x)))
}
