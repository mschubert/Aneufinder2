#' Partition the genome into bins
#'
#' @param variableWidthReference
partitionGenome <- function(seqinfo, binsize=1e6, reads.per.bin=NULL,
                            stepsize=NULL) {
    UseMethod("partitionGenome")
}

partitionGenome.Seqinfo <- function(seqinfo, binsize=1e6, reads.per.bin=NULL,
                                    stepsize=NULL) {
    if (is.null(stepsize))
        stepsize <- binsize
    if (! (binsize >= stepsize && binsize %% stepsize == 0))
        stop("binsize must be a multiple of stepsize")

    # fixed-size bins
    chr.len <- GenomeInfoDb::seqlengths(seqinfo)
    chr.len.floor <- floor(chr.len / stepsize) * stepsize
    # overlapping bins
    if (stepsize < binsize) {
        chr.len.floor <- chr.len.floor - (binsize/stepsize - 1) * stepsize
        chr.len.floor <- pmax(chr.len.floor, 0)
    }

    bins <- GenomicRanges::tileGenome(seqlengths=chr.len.floor, tilewidth=stepsize)
    bins <- unlist(bins, use.names=FALSE)
    if (stepsize < binsize)
        end(bins) <- (start(bins) + binsize - 1)

    GenomeInfoDb::seqinfo(bins) <- seqinfo
    bins
}

#TODO: use custom class here for binned reads(?)
partitionGenome.GRanges <- function(seqinfo, binsize=1e6,
                                    reads.per.bin=NULL, stepsize=NULL) {
    stop("not implemented")
    # seqinfo = GRanges with counts from variable width ref?
    #  would need fixed bins + read counts from readGRanges(variable width ref)

    mediancount <- as.integer(median(bins$counts[bins$counts>0]))
    mc.perstep <- unique(as.integer((1:numsteps) * mediancount / numsteps))

    chrs <- split(bins, GenomicRanges::seqnames(bins))
    margins <- lapply(chrs, function(x) x[seq(mc.perstep, length(x), by=mediancount)])
    margins <- unlist(margins)
#    bins <- 
}

partitionGenome.character <- function(seqinfo, ...) {
    partitionGenome(genome(x), ...)
}

partitionGenome.default <- function(x, ...) {
    stop("Do not know how to create bins from object type ", sQuote(class(x)))
}
