correctGC <- function(reads, method='loess') {
    UseMethod("correctGC")
}

correctGC.GRanges <- function(reads, method='loess') {
    if (method == 'loess') {
        # check: use predict() here instead?
        fit <- stats::loess(counts ~ GC.content, data=as.data.frame(reads))
        correction <- median(counts) / fit$fitted
        reads$counts <- as.integer(round(counts * correction))
        reads$mcounts <- as.integer(round(mcounts * correction))
        reads$pcounts <- as.integer(round(pcounts * correction))
    }

    reads
}

correctGC.list <- function(reads, method='loess') {
    lapply(reads, correctGC, method=method)
}

correctGC.default <- function(reads, method='loess') {
    stop("Do not know how to correct GC for type ", sQuote(class(reads)))
}
