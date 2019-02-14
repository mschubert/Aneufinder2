
# ---------------------------------------------------------------------------------------------------------------------
#                                         variableWidthBins
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 24-09-18
# Last modified: 01-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_variableWidthBins <- function(reads, binned.list){                                                                  # We could also add binsize and stepsize but this should match with binned.list (error prone). Right now we get it from the name of the list item.
  strand(reads)  <- '*'                                                                                                # Sort the reads
  reads          <- sort(reads)
  bins.list      <- list()
  for(ibss in 1:length(binned.list)){                                                                                  # ibss: index bin step size combinations
    bin.and.step <- strsplit(names(binned.list)[ibss],split="binsize_|_stepsize_")                                     # Get binsize and stepsize from name list item.
    binsize      <- as.numeric(bin.and.step[[1]][2])
    stepsize     <- as.numeric(bin.and.step[[1]][3])
    mediancount  <- as.integer(median(binned.list[[ibss]]$counts[binned.list[[ibss]]$counts>0]))                       # Median bincount of fixed width bins.
    numsteps     <- as.integer(binsize / stepsize)                                                                     # In case of overlapping bins (stepsize < binsize).
    mediancount.perstep <- unique(as.integer((1:numsteps) * mediancount / numsteps))
    bins.list.step   <- GRangesList()
    for(istep in 1:length(mediancount.perstep)){                                                                       # Determine consecutive bins for each shift.
      subreads       <- GRangesList()
      skipped.chroms <- character()
      for(chrom in unique(seqnames(reads))){                                                                           # Pick only every mediancount read for each chromosome separately.
        reads.chr    <- reads[seqnames(reads)==chrom]
        if(length(reads.chr) >= mediancount){
          idx        <- seq(mediancount.perstep[istep], length(reads.chr), by=mediancount)
          subreads[[chrom]] <- reads.chr[idx]
        }else{
          skipped.chroms[chrom] <- chrom
        }
      }
      subreads     <- unlist(subreads, use.names=FALSE)
      subreads     <- resize(subreads, width=1)                                                                        # Adjust length of reads to get consecutive bins.
      bins         <- gaps(subreads, start=1L, end=seqlengths(subreads)-1L)                                            # Make new bins. Gaps until seqlengths-1 because we have to add 1 later to get consecutive bins
      bins         <- bins[strand(bins)=='*']
      end(bins)    <- end(bins) + 1
      bins.split   <- split(bins, seqnames(bins))
      if(mediancount.perstep[istep] < mediancount){                                                                    # Remove first bin for smaller than median count steps (overlapping bins).
        bins.split <- endoapply(bins.split, function(x) { x[-1] })
      }
      bins.split   <- endoapply(bins.split, function(x) { x[-length(x)] })                                             # We don't want incomplete bins (number of reads < median count) at the end of each chromosome.
      bins         <- unlist(bins.split, use.names=FALSE)
      bins         <- bins[!seqnames(bins) %in% skipped.chroms]                                                        # Remove skipped chromosomes.
      bins         <- keepSeqlevels(bins, setdiff(seqlevels(bins), skipped.chroms))
      bins.list.step[[as.character(istep)]] <- bins
    }
    bins.list.step <- sort(unlist(bins.list.step, use.names=FALSE))                                                    # Combine steps into one bin object.
    bins.list[[names(binned.list)[ibss]]] <- bins.list.step
  }
  return(bins.list)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------