genome <- function(x, ...) {
    UseMethod("genome")
}

genome.character <- function(x, ...) {
    if (length(x) > 1)
        stop("Can only get genome information for one object at a time")

    if (file.exists(x) && grepl("\\.bam$", x))
        genome(Rsamtools::BamFile(x), ...)
    else {
        lookup <- setNames(c("mm9", "mm10", "hg37", "hg38"),
                           c("GRCm37", "GRCm38", "GRCh37", "GRCh38"))
        adj <- identity
        if (x %in% names(lookup)) {
            adj <- function(x) { GenomeInfoDb::seqlevelsStyle(x) <- "NCBI"; x }
            x <- lookup[[x]]
        }
        adj(genome(GenomeInfoDb::Seqinfo(genome=x), ...))
    }
}

genome.GRanges <- function(x, ...) {
    seqinfo <- GenomeInfoDb::seqinfo(x)
    genome(seqinfo, ...)
}

genome.BamFile <- function(x, ...) {
    seqinfo <- GenomeInfoDb::seqinfo(x)
    genome(seqinfo, ...)
}

genome.Seqinfo <- function(x, chrs=NULL) {
    if (is.null(chrs)) {
        chrs = GenomeInfoDb::standardChromosomes(x)
        chrs = grep("^(chr)?(Y|M(T)?)$", chrs, invert=TRUE, value=TRUE)
    }
    GenomeInfoDb::keepSeqlevels(x, chrs)
}

genome.default <- function(x, ...) {
    stop("Can not parse genome info from type ", sQuote(class(x)))
}
