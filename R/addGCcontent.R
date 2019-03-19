addGCcontent <- function(bins, BSgenome=NULL) {
    UseMethod("addGCcontent")
}

addGCcontent.GRanges <- function(bins, BSgenome=NULL) {
    if (is.null(BSgenome)) {
        assembly <- unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(bins)))
        if (is.na(assembly))
            stop("No genome information in bins and no explicit BSgenome given")
        BSgenome <- grep(assembly, BSgenome::available.genomes(), value=TRUE)[1]
    }
    bs.genome <- getFromNamespace(BSgenome, ns=BSgenome)
    seqlevelsStyle(bs.genome) <- seqlevelsStyle(bins)[1] # [NCBI], Ensembl

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

    bins
}

addGCcontent.list <- function(bins, BSgenome=NULL) {
    lapply(bins, addGCcontent, BSgenome=BSgenome, method=method)
}

addGCcontent.default <- function(bins) {
    stop("Do not know how to get GC for type ", sQuote(class(bins)))
}
