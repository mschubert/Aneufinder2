DNAcopy.findCNVs <- function(binned.data, ID=NULL, CNgrid.start=1.5,
                             strand='*') {

    warlist <- list()
    if (strand=='+')
        select <- 'pcounts'
    else if(strand=='-')
        select <- 'mcounts'
    else if(strand=='*')
        select <- 'counts'

    counts <- mcols(binned.data)[,select]
    if(any(is.na(counts))){
        stop(paste0("ID = ",ID,": NAs found in reads."))
    } else if(all(counts == 0)) {
        wstr = paste0("ID = ",ID,": All counts in data are zero. No DNAcopy done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings     <- warlist
        return(result)
    } else if (any(counts < 0)) {
        wstr = paste0("ID = ",ID,": Some counts in data are negative. No DNAcopy done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings <- warlist
        return(result)
    }

    set.seed(0)
    counts.normal <- (counts+1) / mean(counts+1)
    logcounts <- log2(counts.normal)

    CNA.object <- DNAcopy::CNA(
        genomdat=logcounts,
        maploc=as.numeric(start(binned.data)),
        data.type='logratio',
        # --> We need strand information for the bivariate variant of DNA copy (stacked strands).
        chrom=paste0(as.character(seqnames(binned.data)), "_",
                     as.character(strand(binned.data))))
    CNA.smoothed <- DNAcopy::smooth.CNA(CNA.object)
    CNA.segs <- DNAcopy::segment(CNA.smoothed,verbose=0,min.width=5)
    CNA.segs <- CNA.segs$output
    CNA.segs$chrom <- as.character(CNA.segs$chrom)
    CNA.segs$strand <- substr(CNA.segs$chrom,nchar(CNA.segs$chrom),nchar(CNA.segs$chrom))
    CNA.segs$chrom <- substr(CNA.segs$chrom,1,nchar(CNA.segs$chrom)-2)
    # --> End position is start next - 1
    CNA.segs$loc.end <- c(CNA.segs$loc.start[2:nrow(CNA.segs)] - 1, NA)
    last.ind.chr.strand <- cumsum(rle(paste0(CNA.segs$chrom,"_",CNA.segs$strand))$lengths)
    # --> End position last segment of each chromosome is the length of the chromosome
    CNA.segs$loc.end[last.ind.chr.strand] <- seqlengths(binned.data)[CNA.segs$chrom[last.ind.chr.strand]]

    segs.gr <- GRanges(seqnames = CNA.segs$chrom,
                       ranges = IRanges(start=CNA.segs$loc.start,
                                        end=CNA.segs$loc.end),
                       strand = CNA.segs$strand)
    segs.gr$mean.count <- (2^CNA.segs$seg.mean) * mean(counts+1) - 1
    segs.gr$mean.count[segs.gr$mean.count < 0] <- 0
    # --> Modify bins to contain median count
    seg.ind <- findOverlaps(binned.data,segs.gr,select='first')

    counts.normal <- counts / mean(counts[which(counts > 0)])
    # --> Bit ugly approach!!! And it is definitely not the median!!!
    # --> Why not the median?? I guess you want to get rid of outliers when possible??
    segs.gr$median.count <- sapply(split(counts.normal,seg.ind), function(x) {
        qus <- quantile(x, c(0.01, 0.99))
        y <- x[x >= qus[1] & x <= qus[2]]
        if (sum(y) == 0 | length(y) == 0)
            y <- x
        mean(y)
    })

    counts.median <- segs.gr$median.count[seg.ind]
    CNgrid <- seq(CNgrid.start,6,by=0.01)
    # --> Multiplication of two vectors. First with first, first with second, etc.
    outerRaw <- counts.median %o% CNgrid
    outerDiff <- (outerRaw - round(outerRaw)) ^ 2
    # --> For each muliplication factor the sum of squared differences
    sumOfSquares <- colSums(outerDiff, na.rm=FALSE, dims=1)
    # --> Pick the factor that gave the smallest sum of squared differences.
    # ??????? Is sum of squared differences enough????? For edivisve we also
    # use multiplication factor and the number of segments....!!!!
    CN <- CNgrid[order(sumOfSquares)[1]]
    CN.states <- round(counts.median * CN)
    names(CN.states) <- NULL

    result <- list(ID=ID, bins=binned.data)
    result$bins$state <- factor(paste0(CN.states, '-somy'),
                                levels=paste0(sort(unique(CN.states)),'-somy'))
    result$bins$copy.number <- CN.states
    suppressMessages(
        result$segments <-
            as(collapseBins(as.data.frame(result$bins),
                            column2collapseBy = 'state',
                            columns2drop = 'width',
                            columns2average = c('counts','mcounts','pcounts')),
               'GRanges')
    )
    seqlevels(result$segments) <- seqlevels(result$bins) # --> Correct order from as()
    seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]

    # --> Determine the number of bins for each segment
    chr.strand <- paste0(as.character(seqnames(result$bins)), "_",
                         as.character(strand(result$bins)))
    bin.num <- NULL
    for (chr.st in unique(chr.strand))
        bin.num <- c(bin.num,rle(as.character(result$bins$state[
                which(chr.strand == chr.st)]))$lengths)
    result$segments$num.bins <- bin.num
    mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts",
                                                       "mean.pcounts","state","copy.number")]
    result$weights <- table(result$bins$state) / length(result$bins)
    result$distributions <- assign.distributions(counts=result$bins$counts,
                                                  states=result$bins$state)
    result$warnings <- warlist

    class(result) <- "aneuHMM"
    result
}
