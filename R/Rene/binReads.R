
# ---------------------------------------------------------------------------------------------------------------------
#                                                      binReads
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 24-09-18
# Last modified: 08-01-19
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_binReads <- function(reads, bins){
  reads.plus       <- reads[strand(reads)=='+']
  reads.minus      <- reads[strand(reads)=='-']
  reads.star       <- reads[strand(reads)=='*']
  scounts          <- suppressWarnings(GenomicRanges::countOverlaps(bins,reads.star))
  mcounts          <- suppressWarnings(GenomicRanges::countOverlaps(bins,reads.minus))
  pcounts          <- suppressWarnings(GenomicRanges::countOverlaps(bins,reads.plus))
  counts           <- mcounts + pcounts + scounts
  countmatrix      <- matrix(c(counts,mcounts,pcounts), ncol=3)
  colnames(countmatrix) <- c('counts','mcounts','pcounts')
  mcols(bins)      <- as(countmatrix, 'DataFrame')
  return(bins)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------