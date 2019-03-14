#' Partition the genome into bins
#'
#' @param variableWidthReference
partitionGenome <- function(seqinfo, bin.size=1e6, reads.per.bin=NULL, stepsize=NULL,
                      variable.ref=NULL) {
    UseMethod("partitionGenome")
}

partitionGenome.seqinfo <- function(seqinfo=genome(varable.ref), bin.size=1e6,
                              stepsize=NULL, variable.ref=NULL) {

    if (is.null(variable.ref)) {
        chr.len = GenomeInfoDB::seqlengths(seqinfo)
        if (stepsize == binsize) { # non-overlapping bins
            chr.len.floor <- floor(chr.len / binsize) * binsize
        } else if (stepsize < binsize) { # overlapping bins
            chr.len.floor <- floor(chr.len / stepsize) * stepsize
            chr.len.floor <- chr.len.floor - (((binsize/stepsize)-1)*stepsize)
            chr.len.floor <- pmax(chr.len.floor, 0)
        }

        bins <- unlist(GenomicRanges::tileGenome(seqlengths=chr.len.floor,
                                                 tilewidth=stepsize), use.names=FALSE)

        if (stepsize < binsize)
            end(bins) <- (start(bins) + binsize - 1)
    } else {
        stop("not implemented")
        #TODO: handle variable width ref
    }
}

partitionGenome.character <- function(seqinfo, ...) {
    partitionGenome(genome(x), ...)
}

partitionGenome.default <- function(x, ...) {
    stop("Do not know how to create bins from object type ", sQuote(class(x)))
}
