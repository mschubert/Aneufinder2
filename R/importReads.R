
# ---------------------------------------------------------------------------------------------------------------------
#                                                    importReads
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 01-11-18
# Last modified: 05-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# This R code ...
# ---------------------------------------------------------------------------------------------------------------------

# Patchwork

# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_bam2GRanges <- function(bamfile, bamindex, chrom.lengths, pairedEndReads=FALSE, remove.duplicate.reads=TRUE, min.mapq=10, max.fragment.width=1000, blacklist=NULL){
  ir            <- IRanges(start=rep(1, length(chrom.lengths)), end=chrom.lengths)
  gr            <- GenomicRanges::GRanges(seqnames=names(chrom.lengths), ranges=ir)                                    # Import the file into GRanges
  if(!remove.duplicate.reads){
    sb.param    <- Rsamtools::ScanBamParam(which=range(gr), what='mapq', mapqFilter = min.mapq)
    if(pairedEndReads){
      reads.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=sb.param)
    }else{
      reads.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=sb.param)
    }
  }else{
    sb.param    <- Rsamtools::ScanBamParam(which=range(gr), what='mapq', mapqFilter=min.mapq, flag=scanBamFlag(isDuplicate=FALSE))
    if(pairedEndReads){
      reads.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=sb.param)
    }else{
      reads.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=sb.param)
    }
  }
  if(pairedEndReads){                                                                                                  # Convert to GRanges and filter
    reads <- GenomicAlignments::granges(reads.raw, use.mcols=TRUE, on.discordant.seqnames='drop')                      # treat as one fragment
    reads <- reads[width(reads)<=max.fragment.width]
  }else{
    reads <- GenomicAlignments::granges(reads.raw, use.mcols=TRUE)                                                     # treat as one fragment
  }
  if(!is.null(blacklist)){                                                                                             # Exclude reads falling into blacklisted regions
    overlaps <- findOverlaps(reads, blacklist)
    reads    <- reads[setdiff(1:length(reads), S4Vectors::queryHits(overlaps))]
  }
  if(length(reads) == 0){                                                                                              # RW: Still good to check whether there are still reads left. But maybe at the end??
    if(pairedEndReads){
      stop(paste0("No reads imported. Does your file really contain paired end reads? Try with 'pairedEndReads=FALSE'"))
    }
    stop(paste0('No reads imported! Check your BAM-file ', bamfile))
  }
  reads                 <- keepSeqlevels(reads, as.character(unique(seqnames(reads))))
  na.seqlevels          <- seqlevels(reads)[is.na(seqlengths(reads))]
  reads                 <- reads[seqnames(reads) %in% seqlevels(reads)[!is.na(seqlengths(reads))]]
  reads                 <- keepSeqlevels(reads, as.character(unique(seqnames(reads))))
  if(length(na.seqlevels) > 0){
    warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
  }
  attr(reads, 'ID')     <- basename(bamfile)
  return(reads)
}


RW_bed2GRanges <- function(bedfile, chrom.lengths, remove.duplicate.reads=TRUE, min.mapq=10, blacklist=NULL){
  reads             <- readBed(bedfile, track.line = "auto", remove.unusual = FALSE, zero.based = TRUE)                # Read in bed file
  names(mcols(reads))[1] <- "mapq"
  reads$name        <- NULL
  seqlengths(reads) <- as.numeric(chrom.lengths[names(seqlengths(reads))])                                             # Select chromosomes
  na.seqlevels      <- seqlevels(reads)[is.na(seqlengths(reads))]
  reads             <- reads[which(as.character(seqnames(reads)) %in% names(chrom.lengths))]
  reads             <- keepSeqlevels(reads, as.character(unique(seqnames(reads))))
  if(length(na.seqlevels) > 0){
    warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
  }
  if(!is.na(min.mapq)){
    if(any(is.na(mcols(reads)$mapq))){
      warning(paste0(bedfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NA' to keep all reads."))
      mcols(reads)$mapq[is.na(mcols(reads)$mapq)] <- -1
    }
    reads <- reads[mcols(reads)$mapq >= min.mapq]
  }
  if(remove.duplicate.reads){                                                                                          # Remove duplicate reads
    rows.pos   <- which(as.character(strand(reads)) == '+')                                                            # Changing end position + reads into start position and start position - reads into end position.
    rows.neg   <- which(as.character(strand(reads)) == '-')                                                            # Function duplicated() works uses sequence, strand, start position and width to find duplicates.
    save.start <- start(reads)
    save.end   <- end(reads)
    end(reads)[rows.pos]   <- start(reads)[rows.pos]
    start(reads)[rows.neg] <- end(reads)[rows.neg]
    rows.dup     <- which(duplicated(reads))
    start(reads) <- save.start
    end(reads)   <- save.end
    if(length(rows.dup) > 0){
      reads      <- reads[-rows.dup]
    }
  }
  if(!is.null(blacklist)){                                                                                             # Exclude reads falling into blacklisted regions
    overlaps <- findOverlaps(reads, blacklist)
    reads    <- reads[setdiff(1:length(reads), S4Vectors::queryHits(overlaps))]
  }
  reads                 <- keepSeqlevels(reads, as.character(unique(seqnames(reads))))
  na.seqlevels          <- seqlevels(reads)[is.na(seqlengths(reads))]
  reads                 <- reads[seqnames(reads) %in% seqlevels(reads)[!is.na(seqlengths(reads))]]
  reads                 <- keepSeqlevels(reads, as.character(unique(seqnames(reads))))
  if(length(na.seqlevels) > 0){
    warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
  }
  attr(reads, 'ID')     <- basename(bedfile)
  return(reads)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------