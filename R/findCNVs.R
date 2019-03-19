#' Find copy number variations
#'
#' \code{findCNVs} classifies the binned read counts into several states which
#' represent copy-numbers.
#'
#' @author Aaron Taudt
#' @param method Any combination of \code{c('HMM','dnacopy','edivisive')}.
#' 	 Option \code{method='HMM'} uses a Hidden Markov Model as described in 
#'     doi:10.1186/s13059-016-0971-7 to call copy numbers.
#'   Option \code{'dnacopy'} uses \code{\link[DNAcopy]{segment}} from the 
#'     \pkg{\link[DNAcopy]{DNAcopy}} package to call copy numbers similarly to
#'     the method proposed in doi:10.1038/nmeth.3578, which gives more robust
#'     but less sensitive results compared to the HMM.
#'   Option \code{'edivisive'} (DEFAULT) works like option \code{'dnacopy'} but
#'     uses the \code{\link[ecp]{e.divisive}} function from the \pkg{ecp}
#'     package for segmentation.
#' @inheritParams edivisive.findCNVs
#' @inheritParams HMM.findCNVs
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom stats dgeom dnbinom
#' @export
#'
#' @examples
#' # Get an example BED file with single-cell-sequencing reads
#' bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#' # Bin the data into bin size 1Mp
#' binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                    chromosomes=c(1:19,'X','Y'))
#' # Find copy-numbers
#' model <- findCNVs(binned[[1]])
#' # Check the fit
#' plot(model, type='histogram')
findCNVs <- function(binned, strandseq=FALSE, ID=NULL, method="edivisive",
                     strand='*', R=10, sig.lvl=0.1, eps=0.01, init="standard",
                     max.time=-1, max.iter=1000, num.trials=15,
                     eps.try=max(10*eps, 1), num.threads=1,
                     count.cutoff.quantile=0.999,
                     states=c("zero-inflation",paste0(0:10,"-somy")),
                     most.frequent.state="2-somy",
                     most.frequent.state.strandseq="1-somy", algorithm="EM",
                     initial.params=NULL, verbosity=1) {

    if (!strandseq) {
        if (method == 'HMM') {
            model <- HMM.findCNVs(binned, ID, eps=eps, init=init,
                                  max.time=max.time, max.iter=max.iter,
                                  num.trials=num.trials, eps.try=eps.try,
                                  num.threads=num.threads,
                                  count.cutoff.quantile=count.cutoff.quantile,
                                  strand=strand, states=states,
                                  most.frequent.state=most.frequent.state,
                                  algorithm=algorithm,
                                  initial.params=initial.params,
                                  verbosity=verbosity)
        } else if(method == 'dnacopy') {
            model <- DNAcopy.findCNVs(binned, ID, CNgrid.start=1.5,
                                      strand=strand)
        } else if (method == 'edivisive') {
            model <- edivisive.findCNVs(binned, ID, CNgrid.start=1.5,
                                        strand=strand, R=R, sig.lvl=sig.lvl)
        }
    } else if (strandseq) {
        if (method == 'HMM') {
            model <- biHMM.findCNVs(binned, ID, eps=eps, init=init,
                                    max.time=max.time, max.iter=max.iter,
                                    num.trials=num.trials, eps.try=eps.try,
                                    num.threads=num.threads,
                                    count.cutoff.quantile=count.cutoff.quantile,
                                    states=states,
                                    most.frequent.state.strandseq=most.frequent.state.strandseq,
                                    algorithm=algorithm,
                                    initial.params=initial.params)
        } else if (method == 'dnacopy') {
            model <- biDNAcopy.findCNVs(binned, ID, CNgrid.start=0.5)
        } else if (method == 'edivisive') {
            model <- bi.edivisive.findCNVs(binned, ID, CNgrid.start=0.5, R=R,
                                           sig.lvl=sig.lvl)
        }
    }

#    attr(model, 'call') <- match.call() # this will str(binned)
    model
}
