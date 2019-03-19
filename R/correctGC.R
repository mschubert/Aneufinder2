correctGC <- function(bins, BSgenome=NULL, method='loess') {
    UseMethod("correctGC")
}

correctGC.GRanges <- function(bins, BSgenome=NULL, method='loess') {
    if (is.null(BSgenome)) {
        assembly <- unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(bins)))
        if (is.na(assembly))
            stop("No genome information in bins and no explicit BSgenome given")
        BSgenome <- grep(assembly, BSgenome::available.genomes(), value=TRUE)[1]
    }
    bs.genome <- getFromNamespace(BSgenome, ns=BSgenome)

    # # this needs too much memory: https://support.bioconductor.org/p/89480/
    # views <- Biostrings::Views(bs.genome, bins)
    # bins$GC.content <- Biostrings::letterFrequency(views, letters="GC", as.prob=TRUE)
    chunk2gc <- function(chunk) {
        views <- Biostrings::Views(bs.genome, chunk)
        Biostrings::letterFrequency(views, letters="GC", as.prob=TRUE)
    }
    n.chunks <- 100
    chunks <- split(bins, round(seq(1, n.chunks, length.out=length(bins))))
    bins$GC.content = unlist(lapply(chunks, chunk2gc), use.names=FALSE)

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

correctGC.list <- function(bins, BSgenome=NULL, method='loess') {
    lapply(bins, correctGC, BSgenome=BSgenome, method=method)
}

correctGC.default <- function(bins, method='loess') {
    stop("Do not know how to correct GC for type ", sQuote(class(bins)))
}
