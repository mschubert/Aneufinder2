#' Partition the genome into bins
#'
#' @param variableWidthReference
binGenome <- function(x, bin.size=1e6, reads.per.bin=NULL, stepsize=NULL) {
    UseMethod("binGenome")
}

# fixed bins
binGenome.seqinfo <- function(x, bin.size=1e6, stepsize=NULL) {
    chrom.lengths = GenomeInfoDB::seqlengths(seqinfo)
    if (stepsize == binsize) { # non-overlapping bins
        chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
    } else if (stepsize < binsize) { # overlapping bins
        chrom.lengths.floor <- floor(chrom.lengths / stepsize) * stepsize
        chrom.lengths.floor <- chrom.lengths.floor - (((binsize/stepsize)-1)*stepsize)
        chrom.lengths.floor <- pmax(chrom.lengths.floor, 0)
    }

    bins <- unlist(GenomicRanges::tileGenome(seqlengths=chrom.lengths.floor,
                                             tilewidth=stepsize), use.names=FALSE)

    if (stepsize < binsize)
        end(bins) <- (start(bins) + binsize - 1)

}

# variable bins
binGenome.GRanges <- function(x, seqinfo=genome(x), bin.size=1e6, stepsize=NULL) {
    strand(reads) <- '*'
    reads <- sort(reads)

    mediancount <- as.integer(median(binned.list[[ibss]]$counts[binned.list[[ibss]]$counts>0]))
}

binGenome.character <- function(x, ...) {
    binGenome(readGRanges(x), ...)
}

binGenome.default <- function(x, ...) {
    stop("Do not know how to create bins from object type ", sQuote(class(x)))
}
