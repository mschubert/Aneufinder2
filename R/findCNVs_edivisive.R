#' Find copy number variations (edivisive, univariate)
#'
#' Classify the binned read counts into several states which represent
#' copy-number-variation. The function uses the \code{\link{e.divisive}}
#' function to segment the genome.
#'
#' @param binned.data A \link{GRanges-class} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various
#'   downstream functions. Could be the file name of the \code{binned.data} for
#'   example.
#' @param CNgrid.start Start parameter for the CNgrid variable. Very empiric.
#'   Set to 1.5 for normal data and 0.5 for Strand-seq data.
#' @param strand Find copy-numbers only for the specified strand. One of
#'   \code{c('+', '-', '*')}.
#' @param R method-edivisive: The maximum number of random permutations to use
#'   in each iteration of the permutation test (see \code{\link[ecp]{e.divisive}}).
#'   Increase this value to increase accuracy on the cost of speed.
#' @param sig.lvl method-edivisive: The level at which to sequentially test if
#'   a proposed change point is statistically significant (see
#'   \code{\link[ecp]{e.divisive}}). Increase this value to find more breakpoints.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom ecp e.divisive
edivisive.findCNVs <- function(binned.data, ID=NULL, CNgrid.start=1.5,
                               strand='*', R=10, sig.lvl=0.1) {

    warlist <- list()
    if (strand=='+')
        select <- 'pcounts'
    else if (strand=='-')
        select <- 'mcounts'
    else if (strand=='*')
        select <- 'counts'

    counts <- mcols(binned.data)[,select]
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
        counts.chrom <- counts[chr.rows]
        dim(counts.chrom) <- c(length(counts.chrom),1)
        cp <- ecp::e.divisive(counts.chrom, min.size=2, R=R, sig.lvl=sig.lvl)
        binned.data$cluster[chr.rows] <- cp$cluster + cl
        cl <- cl + length(cp$p.values)
    }

    # --> Modify bins to contain mean count
    counts.normal <- counts / mean(counts[which(counts > 0)])
    cnmean <- sapply(split(counts.normal,binned.data$cluster), function(x) {
        qus <- quantile(x, c(0.01, 0.99))
        y <- x[x >= qus[1] & x <= qus[2]]
        if(sum(y) == 0 | length(y) == 0)
            y <- x
        mean(y)
    })
    counts.normal.mean <- cnmean[as.character(binned.data$cluster)]
    CNgrid <- seq(CNgrid.start,6,by=0.01) # --> Determine copy number
    outerRaw <- counts.normal.mean %o% CNgrid
    outerDiff <- (outerRaw - round(outerRaw)) ^ 2
    sumOfSquares <- colSums(outerDiff,na.rm=FALSE,dims=1)
    CN <- CNgrid[order(sumOfSquares)][1]
    CN.states <- round(counts.normal.mean * CN)
    names(CN.states) <- NULL

    result <- list(ID=ID, bins=binned.data)
    result$bins$state <- factor(paste0(CN.states,'-somy'),
                                levels=paste0(sort(unique(CN.states)),'-somy'))
    result$bins$copy.number <- CN.states
    suppressMessages(result$segments <-
        as(collapseBins(as.data.frame(result$bins),
                        column2collapseBy = 'state',
                        columns2drop = 'width',
                        columns2average = c('counts','mcounts','pcounts')),
        'GRanges')
    )
    seqlevels(result$segments) <- seqlevels(result$bins) # --> Correct order from as()
    seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]

    bin.num <- NULL # --> Determine the number of bins for each segment
    for (chr in unique(as.character(seqnames(result$bins))))
        bin.num <- c(bin.num, rle(as.character(result$bins$state[
                which(as.character(seqnames(result$bins)) == chr)]))$lengths)
    result$segments$num.bins <- bin.num
    mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts",
                                                     "mean.pcounts","state","copy.number")]
    result$bins$cluster <- NULL
    result$weights <- table(result$bins$state) / length(result$bins)
    result$distributions <- assign.distributions(counts=result$bins$counts,
                                                 states=result$bins$state)
    result$warnings <- warlist

    class(result) <- "aneuHMM"
    result
}
