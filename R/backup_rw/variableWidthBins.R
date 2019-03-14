variableWidthBins <- function(reads, binned.list) {
    # We could also add binsize and stepsize but this should match with
    # binned.list (error prone). Right now we get it from the name of the list
    # item.
    strand(reads) <- '*' # Sort the reads
    reads <- sort(reads)
    bins.list <- list()

    # ibss: index bin step size combinations
    for (ibss in 1:length(binned.list)) {
        # Get binsize and stepsize from name list item.
        bin.and.step <- strsplit(names(binned.list)[ibss],split="binsize_|_stepsize_")

        binsize <- as.numeric(bin.and.step[[1]][2])
        stepsize <- as.numeric(bin.and.step[[1]][3])

        # Median bincount of fixed width bins.
        mediancount <- as.integer(median(binned.list[[ibss]]$counts[binned.list[[ibss]]$counts>0]))
        # In case of overlapping bins (stepsize < binsize).
        numsteps <- as.integer(binsize / stepsize)
        mediancount.perstep <- unique(as.integer((1:numsteps) * mediancount / numsteps))
        bins.list.step <- GRangesList()

        # Determine consecutive bins for each shift.
        for (istep in 1:length(mediancount.perstep)) {
            subreads <- GRangesList()
            skipped.chroms <- character()

            # Pick only every mediancount read for each chromosome separately.
            for (chrom in unique(seqnames(reads))) {
                reads.chr <- reads[seqnames(reads)==chrom]
                if (length(reads.chr) >= mediancount){
                    idx <- seq(mediancount.perstep[istep], length(reads.chr), by=mediancount)
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
            if (mediancount.perstep[istep] < mediancount) {
                bins.split <- endoapply(bins.split, function(x) { x[-1] })
            }

            # We don't want incomplete bins (number of reads < median count) at the
            # end of each chromosome.
            bins.split <- endoapply(bins.split, function(x) { x[-length(x)] })
            bins <- unlist(bins.split, use.names=FALSE)
            bins <- bins[!seqnames(bins) %in% skipped.chroms]
            bins <- keepSeqlevels(bins, setdiff(seqlevels(bins), skipped.chroms))
            bins.list.step[[as.character(istep)]] <- bins
        }

        # Combine steps into one bin object.
        bins.list.step <- sort(unlist(bins.list.step, use.names=FALSE))
        bins.list[[names(binned.list)[ibss]]] <- bins.list.step
    }

    bins.list
}
