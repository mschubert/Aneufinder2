getGCContentBins <- function(bins.list, GC.BSgenome){
    bins.list.gc <- list()
    # ibss: index bin step size combinations
    for (ibss in 1:length(bins.list)) {
        bin <- bins.list[[ibss]]
        GC.content <- list()

        for (chr in unique(as.character(seqnames(bin)))) {
            view <- Biostrings::Views(GC.BSgenome[[chr]], ranges(bin)[seqnames(bin)==chr])
            freq <- Biostrings::alphabetFrequency(view, as.prob = TRUE, baseOnly=TRUE)
            GC.content[[as.character(chr)]] <- rowSums(freq[, c("G","C"), drop=FALSE])
        }

        GC.content <- unlist(GC.content)
        bin$GC <- GC.content
        bins.list.gc[[names(bins.list)[ibss]]] <- bin
    }

    bins.list.gc
}

correctGC <- function(binned.gc, method='loess') {
    counts <- binned.gc$counts
    mcounts <- binned.gc$mcounts
    pcounts <- binned.gc$pcounts
    GC <- binned.gc$GC

    if (method == 'loess') {
        mean.counts.global <- median(counts)
        fit <- stats::loess(counts ~ GC)
        correction.factor <- mean.counts.global / fit$fitted
        counts <- counts * correction.factor
        mcounts <- mcounts * correction.factor
        pcounts <- pcounts * correction.factor
    }

    binned.gc$counts <- as.integer(round(counts))
    binned.gc$mcounts <- as.integer(round(mcounts))
    binned.gc$pcounts <- as.integer(round(pcounts))
    binned.gc
}
