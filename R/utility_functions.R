transCoord <- function(gr) {
  cum.seqlengths        <- cumsum(as.numeric(seqlengths(gr)))
  cum.seqlengths.0      <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- seqlevels(gr)
  gr$start.genome       <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  gr$end.genome         <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  return(gr)
}

