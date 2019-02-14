
# ---------------------------------------------------------------------------------------------------------------------
#                                         FixedWidthBins
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 21-09-18
# Last modified: 01-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_fixedWidthBins <- function(chrom.lengths=NULL, binsizes=1e6, stepsizes=NULL){                                       # Could also add a list with chromosomes that we selected. But not necessary if chromosome lengths only contains these selected ones.
  bins.list  <- list()
  for(ibss in 1:length(binsizes)){                                                                                     # ibss: index bin step size combinations
    binsize  <- binsizes[ibss]
    stepsize <- stepsizes[ibss]
    if(stepsize == binsize){                                                                                           # Non-overlapping bins.
      chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
    }else if(stepsize < binsize){                                                                                      # Overlapping bins.
      chrom.lengths.floor <- floor(chrom.lengths / stepsize) * stepsize
      chrom.lengths.floor <- chrom.lengths.floor - (((binsize/stepsize)-1)*stepsize)
      chrom.lengths.floor[which(chrom.lengths.floor < 0)] <- 0                                                         # The last step can result in negative chromosome lengths.
    }
    bins <- unlist(GenomicRanges::tileGenome(seqlengths=chrom.lengths.floor, tilewidth=stepsize), use.names=FALSE)
    seqlengths(bins) <- chrom.lengths
    if(stepsize < binsize){                                                                                            # Change binsize to true binsize (for creating overlapping bins).
      end(bins) <- (start(bins) + binsize - 1)
    }
    if (any(width(bins) != binsize)){                                                                                  # None of the bins have the correct size.
      stop("tileGenome failed")
    }
    bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE), '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins
    skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
    if(length(skipped.chroms)>0) {
      warning("Chromosomes smaller than binsize ", binsize, " (skipped): ", paste0(skipped.chroms, collapse=', '))
    }
  }
  return(bins.list)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------