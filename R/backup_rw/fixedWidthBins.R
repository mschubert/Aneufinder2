fixedWidthBins <- function(chrom.lengths=NULL, binsizes=1e6, stepsizes=NULL) {
    # Could also add a list with chromosomes that we selected. But not
    # necessary if chromosome lengths only contains these selected ones.
    bins.list  <- list()

    # ibss: index bin step size combinations
    for (ibss in 1:length(binsizes)) {
        binsize  <- binsizes[ibss]
        stepsize <- stepsizes[ibss]

        if(stepsize == binsize) {
            # Non-overlapping bins
            chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
        } else if(stepsize < binsize) {
            # overlapping bins
            chrom.lengths.floor <- floor(chrom.lengths / stepsize) * stepsize
            chrom.lengths.floor <- chrom.lengths.floor - (((binsize/stepsize)-1)*stepsize)
            # The last step can result in negative chromosome lengths.
            chrom.lengths.floor[which(chrom.lengths.floor < 0)] <- 0
        }

        bins <- unlist(GenomicRanges::tileGenome(seqlengths=chrom.lengths.floor,
                                                 tilewidth=stepsize), use.names=FALSE)
        seqlengths(bins) <- chrom.lengths

        # Change binsize to true binsize (for creating overlapping bins).
        if (stepsize < binsize)
            end(bins) <- (start(bins) + binsize - 1)

        # None of the bins have the correct size.
        if (any(width(bins) != binsize))
            stop("tileGenome failed")

        bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE),
                          '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins
        skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))

        if (length(skipped.chroms) > 0)
            warning("Chromosomes smaller than binsize ", binsize, " (skipped): ",
                    paste0(skipped.chroms, collapse=', '))
    }

    bins.list
}
