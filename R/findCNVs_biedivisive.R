#' Find copy number variations (edivisive, bivariate)
#'
#' Classify the binned read counts into several states which represent
#' copy-number-variation. The function uses the \code{\link{e.divisive}}
#' function to segment the genome.
#'
#' @inheritParams edivisive.findCNVs
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom ecp e.divisive
bi.edivisive.findCNVs <- function(binned.data, ID=NULL, CNgrid.start=0.5, R=10,
                                  sig.lvl=0.1) {

    warlist <- list()
    counts <- as.matrix(mcols(binned.data)[,c('mcounts','pcounts')])

    # --> Check if there are counts in the data  --> Check with other bi functions!!!
    if (any(is.na(counts))) {
        stop(paste0("ID = ", ID, ": NAs found in reads."))
    } else if(all(counts == 0)) {
        wstr = paste0("ID = ", ID, ": All counts in data are zero. No Edivisive done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings <- warlist
        return(result)
    } else if (any(counts < 0)) {
        wstr = paste0("ID = ", ID, ": Some counts in data are negative. No Edivisive done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings <- warlist
        return(result)
    }

    set.seed(0)
    binned.data$cluster <- NA
    cl <- 0
    for (chrom in unique(as.character(seqnames(binned.data)))) {
        chr.rows <- which(as.character(seqnames(binned.data)) == chrom)
        counts.chrom <- counts[chr.rows,]
        cp <- ecp::e.divisive(counts.chrom, min.size=2, R=R, sig.lvl=sig.lvl)
        binned.data$cluster[chr.rows] <- cp$cluster + cl
        cl <- cl + length(cp$p.values)
    }

    # --> Modify bins to contain median count
    counts.normal <- counts / mean(counts[which(counts > 0)])
    cnmean.m <- numeric()
    cnmean.p <- numeric()
    for (i1 in 1:max(binned.data$cluster)) {
        x                    <- counts.normal[binned.data$cluster==i1,,drop=FALSE]
        qus                  <- quantile(x, c(0.01, 0.99))
        # --> The counts of both strands need to be within the quantile range
        within.quantile      <- apply(x,2,function(z){z >= qus[1] & z <= qus[2]})
        dim(within.quantile) <- dim(x) # --> Keep the dimensions of thw matrix
        dimnames(within.quantile) <- dimnames(x)
        within.quantile      <- within.quantile[,'mcounts'] & within.quantile[,'pcounts']
        y                    <- x[within.quantile,,drop=FALSE]
        if (sum(y) == 0 | length(y)==0)
            y <- x
        mu <- colMeans(y)
        cnmean.m[as.character(i1)] <- mu['mcounts']
        cnmean.p[as.character(i1)] <- mu['pcounts']
    }

    counts.normal.mean.m  <- cnmean.m[as.character(binned.data$cluster)]
    counts.normal.mean.p  <- cnmean.p[as.character(binned.data$cluster)]
    counts.normal.mean.stacked <- c(counts.normal.mean.m,counts.normal.mean.p)
    CNgrid                <- seq(CNgrid.start,6,by=0.01)
    outerRaw              <- counts.normal.mean.stacked %o% CNgrid
    outerDiff             <- (outerRaw - round(outerRaw)) ^ 2
    sumOfSquares          <- colSums(outerDiff,na.rm=FALSE,dims=1)
    CN                    <- CNgrid[order(sumOfSquares)][1]
    CN.states             <- round(counts.normal.mean.stacked * CN)
    names(CN.states)      <- NULL

    result                <- list(ID=ID, bins=binned.data)
    result$bins$mcopy.number <- CN.states[1:(length(CN.states)/2)]
    result$bins$pcopy.number <- CN.states[((length(CN.states)/2)+1):length(CN.states)]
    result$bins$copy.number  <- result$bins$mcopy.number + result$bins$pcopy.number
    lev                   <- paste0(0:max(result$bins$copy.number),'-somy')
    result$bins$state     <- factor(paste0(result$bins$copy.number,'-somy'),levels=lev)
    result$bins$mstate    <- factor(paste0(result$bins$mcopy.number,'-somy'),levels=lev)
    result$bins$pstate    <- factor(paste0(result$bins$pcopy.number,'-somy'),levels=lev)
    result$bins$combi     <- paste(result$bins$mcopy.number,result$bins$pcopy.number)
    suppressMessages(result$segments <-
        as(collapseBins(as.data.frame(result$bins), column2collapseBy='combi',
                        columns2drop='width',
                        columns2average=c('counts','mcounts','pcounts')),
           'GRanges')
    )
    seqlevels(result$segments) <- seqlevels(result$bins) # --> Correct order from as()
    seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]

    bin.num <- NULL # --> Determine the number of bins for each segment
    for (chr in unique(as.character(seqnames(result$bins))))
        bin.num <- c(bin.num, rle(result$bins$combi[
                which(as.character(seqnames(result$bins)) == chr)])$lengths)
    result$segments$num.bins <- bin.num
    mcols(result$segments) <-
        mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts","mean.pcounts","state",
                                 "mstate","pstate","copy.number","mcopy.number","pcopy.number")]
    mcols(result$bins) <- mcols(result$bins)[c("counts","mcounts","pcounts","GC","state","mstate",
                                               "pstate","copy.number","mcopy.number","pcopy.number")]
    result$weights <- table(result$bins$state) / length(result$bins)
    dist.strand <- assign.distributions(counts=c(result$bins$mcounts,result$bins$pcounts),
         states=c(as.character(result$bins$mstate),as.character(result$bins$pstate)))
    result$distributions <- list(minus = dist.strand, plus = dist.strand)
    result$distributions$both <- assign.distributions(counts=result$bins$counts,states=result$bins$state)
    result$univariateParams <- list(weights=table(c(as.character(result$bins$mstate),
                               as.character(result$bins$pstate))) / (length(result$bins)*2))
    result$warnings <- warlist

    class(result) <- "aneuBiHMM"
    result
}
