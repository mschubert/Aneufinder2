findCNVs <- function(strandseq=FALSE, binned, ID=NULL, method="edivisive",
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

    attr(model, 'call') <- match.call()
    model
}
