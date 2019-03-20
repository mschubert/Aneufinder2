adjustPartitions <- function(bins, reads) {
    bins <- binReads(reads, bins)
    strand(reads) <- '*'
    reads <- sort(reads)

    binsize <- attr(bins, 'binsize')
    stepsize <- attr(bins, 'stepsize')
    numsteps <- binsize / stepsize

    mediancount <- as.integer(median(bins$counts[bins$counts>0]))
    mc.perstep <- unique(as.integer((1:numsteps) * mediancount / numsteps))

    # Determine consecutive bins for each shift.
    for (istep in 1:length(mc.perstep)) {
        subreads <- GRangesList()
        skipped.chroms <- character()

        # Pick only every mediancount read for each chromosome separately.
        for (chrom in unique(seqnames(reads))) {
            reads.chr <- reads[seqnames(reads)==chrom]
            if (length(reads.chr) >= mediancount){
                idx <- seq(mc.perstep[istep], length(reads.chr), by=mediancount)
                subreads[[chrom]] <- reads.chr[idx]
            } else {
                skipped.chroms[chrom] <- chrom
            }
        }

        subreads <- unlist(subreads, use.names=FALSE)
        # Adjust length of reads to get consecutive bins.
        subreads <- resize(subreads, width=1)
        # Make new bins. Gaps until seqlengths-1 because we have to add 1 later
        # to get consecutive bins
        bins <- gaps(subreads, start=1L, end=seqlengths(subreads)-1L)
        bins <- bins[strand(bins)=='*']
        end(bins) <- end(bins) + 1
        bins.split <- split(bins, seqnames(bins))

        # Remove first bin for smaller than median count steps (overlapping bins).
        if (mc.perstep[istep] < mediancount) {
            bins.split <- endoapply(bins.split, function(x) { x[-1] })
        }

        # We don't want incomplete bins (number of reads < median count) at the
        # end of each chromosome.
        bins.split <- endoapply(bins.split, function(x) { x[-length(x)] })
        bins <- unlist(bins.split, use.names=FALSE)
        bins <- bins[!seqnames(bins) %in% skipped.chroms]
        bins <- keepSeqlevels(bins, setdiff(seqlevels(bins), skipped.chroms))
    }

    bins
}
