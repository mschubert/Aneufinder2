
# ---------------------------------------------------------------------------------------------------------------------
#                                                      correctGC
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 26-09-18
# Last modified: 02-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_getGCContentBins <- function(bins.list,GC.BSgenome){
  bins.list.gc          <- list()
  for(ibss in 1:length(bins.list)){                                                                                    # ibss: index bin step size combinations
    bin                 <- bins.list[[ibss]]
    GC.content          <- list()
    for (chr in unique(as.character(seqnames(bin)))){
      view              <- Biostrings::Views(GC.BSgenome[[chr]], ranges(bin)[seqnames(bin)==chr])
      freq              <- Biostrings::alphabetFrequency(view, as.prob = TRUE, baseOnly=TRUE)
      GC.content[[as.character(chr)]] <- rowSums(freq[, c("G","C"), drop=FALSE])
    }
    GC.content          <- unlist(GC.content)
    bin$GC              <- GC.content
    bins.list.gc[[names(bins.list)[ibss]]] <- bin
  }
  return(bins.list.gc)
}


RW_correctGC <- function(binned.gc,method='loess'){
  counts                <- binned.gc$counts
  mcounts               <- binned.gc$mcounts
  pcounts               <- binned.gc$pcounts
  GC                    <- binned.gc$GC
  if(method == 'loess'){
    mean.counts.global  <- median(counts)
    fit                 <- stats::loess(counts ~ GC)
    correction.factor   <- mean.counts.global / fit$fitted
    counts              <- counts * correction.factor
    mcounts             <- mcounts * correction.factor
    pcounts             <- pcounts * correction.factor
  }
  binned.gc$counts      <- as.integer(round(counts))
  binned.gc$mcounts     <- as.integer(round(mcounts))
  binned.gc$pcounts     <- as.integer(round(pcounts))
  return(binned.gc)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------