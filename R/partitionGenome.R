#' Partition the genome into bins
#'
#' @param variableWidthReference
partitionGenome <- function(seqinfo, bin.size=1e6, reads.per.bin=NULL,
                            stepsize=NULL) {
    UseMethod("partitionGenome")
}

partitionGenome.Seqinfo <- function(seqinfo, bin.size=1e6, stepsize=NULL) {
    if (is.null(stepsize))
        stepsize <- bin.size
    if (! (bin.size >= stepsize && bin.size %% stepsize == 0))
        stop("bin.size must be a multiple of stepsize")

    # fixed-size bins
    chr.len <- GenomeInfoDb::seqlengths(seqinfo)
    chr.len.floor <- floor(chr.len / stepsize) * stepsize
    # overlapping bins
    if (stepsize < bin.size) {
        chr.len.floor <- chr.len.floor - (bin.size/stepsize - 1) * stepsize
        chr.len.floor <- pmax(chr.len.floor, 0)
    }

    bins <- GenomicRanges::tileGenome(seqlengths=chr.len.floor, tilewidth=stepsize)
    bins <- unlist(bins, use.names=FALSE)
    if (stepsize < bin.size)
        end(bins) <- (start(bins) + bin.size - 1)

    GenomeInfoDb::seqinfo(bins) <- seqinfo
    bins
}

partitionGenome.GRanges <- function(seqinfo, bin.size=1e6, stepsize=NULL) {
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
