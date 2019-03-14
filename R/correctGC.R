correctGC <- function(bins, method='loess') {
    UseMethod("correctGC")
}

correctGC.GRanges <- function(bins, method='loess') {
    assembly <- unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(bins)))
    bs.name <- grep(assembly, BSgenome::available.genomes(), value=TRUE)[1]
    bs.genome <- getFromNamespace(bs.name, ns=bs.name)

    # check memory here; cf: https://support.bioconductor.org/p/89480/
    views <- Biostrings::Views(bs.genome, bins)
    bins$GC.content <- Biostrings::letterFrequency(views, letters="GC", as.prob=TRUE)

    if (method == 'loess') {
        # check: use predict() here instead?
        fit <- stats::loess(counts ~ GC.content, data=as.data.frame(bins))
        correction <- median(counts) / fit$fitted
        bins$counts <- as.integer(round(counts * correction))
        bins$mcounts <- as.integer(round(mcounts * correction))
        bins$pcounts <- as.integer(round(pcounts * correction))
    }

    bins
}

correctGC.list <- function(bins, method='loess') {
    lapply(bins, correctGC, method=method)
}

correctGC.default <- function(bins, method='loess') {
    stop("Do not know how to correct GC for type ", sQuote(class(bins)))
}
