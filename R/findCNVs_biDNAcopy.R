#' Find copy number variations (DNAcopy, bivariate)
#'
#' \code{biDNAcopy.findCNVs} classifies the binned read counts into several
#' states which represent copy-number-variation using read count information
#' from both strands.
#'
#' @param binned.data A \link{GRanges-class} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various
#'   downstream functions. Could be the file name of the \code{binned.data} for
#'   example.
#' @param CNgrid.start Start parameter for the CNgrid variable. Very empiric.
#'   Set to 1.5 for normal data and 0.5 for Strand-seq data.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom DNAcopy CNA smooth.CNA
biDNAcopy.findCNVs <- function(binned.data,ID=NULL, CNgrid.start=0.5) {
    warlist <- list()
    counts <- matrix(c(mcols(binned.data)[,'mcounts'],
                       mcols(binned.data)[,'pcounts']), ncol=2,
                     dimnames=list(bin=1:length(binned.data),
                                   strand=c('minus','plus')))

    if (any(is.na(counts))) {
        stop(paste0("ID = ",ID,": NAs found in reads."))
    } else if(all(counts == 0)) {
        wstr = paste0("ID = ",ID,": All counts in data are zero. No DNAcopy done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings     <- warlist
        return(result)
    } else if (any(counts < 0)) {
        wstr = paste0("ID = ",ID,": Some counts in data are negative. No DNAcopy done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings     <- warlist
        return(result)
    }

    # STACK THE STRANDS AND RUN DNACOPY
    binned.data.minus     <- binned.data
    strand(binned.data.minus) <- '-'
    binned.data.minus$counts <- binned.data.minus$mcounts
    binned.data.plus      <- binned.data
    strand(binned.data.plus) <- '+'
    binned.data.plus$counts <- binned.data.plus$pcounts
    binned.data.stacked   <- c(binned.data.minus,binned.data.plus)
    message("Running DNAcopy")
    model.stacked         <- DNAcopy.findCNVs(binned.data.stacked,ID,CNgrid.start=CNgrid.start)

    result                <- list(ID=ID, bins=binned.data)
    result$bins$state     <- NA
    result$bins$mstate    <- model.stacked$bins$state[as.logical(model.stacked$bins@strand=='-')]
    result$bins$pstate    <- model.stacked$bins$state[as.logical(model.stacked$bins@strand=='+')]
    result$bins$copy.number <- NA
    result$bins$mcopy.number <- model.stacked$bins$copy.number[as.logical(model.stacked$bins@strand=='-')]
    result$bins$pcopy.number <- model.stacked$bins$copy.number[as.logical(model.stacked$bins@strand=='+')]
    result$bins$copy.number <- result$bins$mcopy.number + result$bins$pcopy.number
    bs.state              <- paste0(result$bins$copy.number,"-somy")
    result$bins$state     <- factor(bs.state, levels=sort(unique(c(as.character(result$bins$mstate),
                             as.character(result$bins$pstate),bs.state))))
    result$bins$combi     <- paste(result$bins$mcopy.number,result$bins$pcopy.number)
    suppressMessages(result$segments <-
        as(collapseBins(as.data.frame(result$bins), column2collapseBy='combi',
                        columns2drop='width',
                        columns2average=c('counts','mcounts','pcounts')),
           'GRanges')
    )
    seqlevels(result$segments) <- seqlevels(result$bins) # --> Correct order from as()
    seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]

    bin.num <- NULL
    for(chr in unique(as.character(seqnames(result$bins))))
        bin.num <- c(bin.num,rle(result$bins$combi[
            which(as.character(seqnames(result$bins)) == chr)])$lengths)
    result$segments$num.bins <- bin.num
    mcols(result$segments) <-
        mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts","mean.pcounts","state",
                                 "mstate","pstate","copy.number","mcopy.number","pcopy.number")]
    result$bins$combi     <- NULL
    result$weights        <- table(result$bins$state) / length(result$bins)
    result$distributions  <- list(minus = model.stacked$distributions, plus = model.stacked$distributions)
    result$distributions$both <- assign.distributions(counts=result$bins$counts, states=result$bins$state)
    result$univariateParams <- list(weights=model.stacked$weights)
    result$warnings       <- model.stacked$warnings

    class(result) <- "aneuBiHMM"
    result
}
