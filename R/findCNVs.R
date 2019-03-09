
# ---------------------------------------------------------------------------------------------------------------------
#                                                      findCNVs
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 09-11-18
# Last modified: 23-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_findCNVs <- function(strandseq=FALSE, binned, ID=NULL, method="edivisive", strand='*', R=10, sig.lvl=0.1, eps=0.01,
                 init="standard", max.time=-1, max.iter=1000, num.trials=15, eps.try=max(10*eps, 1), num.threads=1,
                 count.cutoff.quantile=0.999, states=c("zero-inflation",paste0(0:10,"-somy")),
                 most.frequent.state="2-somy", most.frequent.state.strandseq="1-somy", algorithm="EM",
                 initial.params=NULL, verbosity=1){
  if(!strandseq){
    if(method == 'HMM'){
      model <- RW_HMM.findCNVs(binned, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter,
                 num.trials=num.trials, eps.try=eps.try, num.threads=num.threads,
                 count.cutoff.quantile=count.cutoff.quantile, strand=strand, states=states,
                 most.frequent.state=most.frequent.state, algorithm=algorithm, initial.params=initial.params,
                 verbosity=verbosity)
    }else if(method == 'dnacopy'){
      model <- RW_DNAcopy.findCNVs(binned, ID, CNgrid.start=1.5, strand=strand)
    }else if (method == 'edivisive'){
      model <- RW_edivisive.findCNVs(binned, ID, CNgrid.start=1.5, strand=strand, R=R, sig.lvl=sig.lvl)
    }
  }else if(strandseq){
    if(method == 'HMM'){
      model <- RW_biHMM.findCNVs(binned, ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter,
                 num.trials=num.trials, eps.try=eps.try, num.threads=num.threads,
                 count.cutoff.quantile=count.cutoff.quantile, states=states,
                 most.frequent.state.strandseq=most.frequent.state.strandseq, algorithm=algorithm,
                 initial.params=initial.params)
    }else if (method == 'dnacopy'){
      model <- RW_biDNAcopy.findCNVs(binned, ID, CNgrid.start=0.5)
    }else if (method == 'edivisive'){
      model <- RW_bi.edivisive.findCNVs(binned, ID, CNgrid.start=0.5, R=R, sig.lvl=sig.lvl)
    }
  }
  attr(model, 'call') <- match.call()
  return(model)
}


RW_HMM.findCNVs <- function(binned.data, ID=NULL, eps=0.01, init="standard", max.time=-1, max.iter=-1, num.trials=1,   # ==============================================================================
                     eps.try=NULL, num.threads=1, count.cutoff.quantile=0.999, strand='*',                             # GET THE COUNTS AND CHECK THE DATA
                     states=c("zero-inflation",paste0(0:10,"-somy")), most.frequent.state="2-somy", algorithm="EM",    # ==============================================================================
                     initial.params=NULL, verbosity=1){
  on.exit(.C("C_univariate_cleanup", PACKAGE = 'AneuFinder'))                                                          # --> Define cleanup behaviour
  warlist               <- list()
  if(!is.null(initial.params)){
    init                <- 'initial.params'
  }
  if(num.trials == 1){
    eps.try             <- eps
  }
  if(strand=='+'){                                                                                                     # --> Get counts
    select              <- 'pcounts'
  }else if (strand=='-'){
    select              <- 'mcounts'
  }else if (strand=='*'){
    select              <- 'counts'
  }
  counts                <- mcols(binned.data)[,select]
  count.cutoff          <- ceiling(quantile(counts, count.cutoff.quantile))                                            # --> Filter high counts out, makes HMM faster
  rows.higher           <- which(counts > count.cutoff)
  counts[rows.higher]   <- count.cutoff
  if(length(rows.higher) > 0){
    message(paste0("Replaced read counts > ",count.cutoff," (quantile ",count.cutoff.quantile,") by ",count.cutoff,
	" in ",length(rows.higher)," bins. Set option 'count.cutoff.quantile=1' to disable this filtering.",
    " This filtering was done to enhance performance."))
  }
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data, otherwise HMM will blow up
    stop(paste0("ID = ",ID,": NAs found in reads."))
  }else if(!any(counts != 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No HMM done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No HMM done."))
    result$warnings     <- warlist
    return(result)
  }
  
  numbins               <- length(binned.data)
  numstates             <- length(states)
  
  inistates             <- initializeStates(states)                                                                    # --> Assign variables
  state.labels          <- inistates$states
  state.distributions   <- inistates$distributions
  multiplicity          <- inistates$multiplicity
  dependent.states.mask <- (state.labels != 'zero-inflation') & (state.labels != '0-somy')
  
  algorithm             <- factor(algorithm,levels=c('baumWelch','viterbi','EM'))                                      # --> This factor is later on changed into a number.

  modellist <- list()                                                                                                  # ==============================================================================
  for(i_try in 1:num.trials){                                                                                          # RUN MODEL WITH EACH TIME DIFFERENT INITIAL PARAMETERS
    if (verbosity >= 1){message(paste0("Trial ",i_try," / ",num.trials))}                                              # ==============================================================================
    if(init == 'initial.params'){                                                                                      # --> Initial parameters
      A.initial         <- initial.params$transitionProbs
      proba.initial     <- initial.params$startProbs
      size.initial      <- initial.params$distributions[,'size']
      prob.initial      <- initial.params$distributions[,'prob']
      size.initial[is.na(size.initial)] <- 0
      prob.initial[is.na(prob.initial)] <- 0
    }else if(init == 'random'){
      A.initial         <- matrix(stats::runif(numstates^2),ncol=numstates)
      A.initial         <- sweep(A.initial,1,rowSums(A.initial),"/")			
      proba.initial     <- stats::runif(numstates)
      size.initial      <- stats::runif(1, min=0, max=100) * cumsum(dependent.states.mask)                             # --> Distributions for dependent states
      prob.initial      <- stats::runif(1) * dependent.states.mask                                                     # --> Do we assume a mean monosomy between 0 and 100? Isn't this a bit low????
      index             <- which(state.labels == '0-somy')                                                             # --> Assign initials for the 0-somy distribution
      size.initial[index] <- 1
      prob.initial[index] <- 0.5
    }else if(init == 'standard'){
      A.initial         <- matrix((0.1/(numstates-1)), ncol=numstates, nrow=numstates)
      for(i in 1:numstates){
        A.initial[i,i]  <- 0.9
      }
      proba.initial     <- rep(1/numstates, numstates)
      max.counts        <- as.integer(names(which.max(table(counts[counts>0]))))                                       # --> Set initial mean of most.frequent.state distribution to max of count histogram
      divf              <- max(multiplicity[most.frequent.state], 1)
      mean.initial.monosomy <- max.counts/divf
      if(is.na(mean.initial.monosomy)) {
        mean.initial.monosomy <- 1
      }
      mean.initial      <- mean.initial.monosomy * cumsum(dependent.states.mask)                                       # --> Do you not assume here that the states are consecutive (e.g. 1,2,3 not 1,3,4,...)?
      var.initial       <- mean.initial*2
      size.initial      <- rep(0,numstates)
      prob.initial      <- rep(0,numstates)
      mask              <- dependent.states.mask
      size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
      prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
      index             <- which(state.labels == '0-somy')                                                             # --> Assign initials for the 0-somy distribution
      size.initial[index] <- 1
      prob.initial[index] <- 0.5
    }
    hmm                 <- .C("C_univariate_hmm",
                              counts = as.integer(counts),                                                             # --> int* O
                              num.bins = as.integer(numbins),                                                          # --> int* T
                              num.states = as.integer(numstates),                                                      # --> int* N
                              state.labels = as.integer(state.labels),                                                 # --> int* state_labels
                              size = double(length=numstates),                                                         # --> double* size
                              prob = double(length=numstates),                                                         # --> double* prob
                              num.iterations = as.integer(max.iter),                                                   # --> int* maxiter
                              time.sec = as.integer(max.time),                                                         # --> double* maxtime
                              loglik.delta = as.double(eps.try),                                                       # --> double* eps
                              maxPosterior = double(length=numbins),                                                   # --> double* maxPosterior
                              states = integer(length=numbins),                                                        # --> int* states
                              A = double(length=numstates*numstates),                                                  # --> double* A
                              proba = double(length=numstates),                                                        # --> double* proba
                              loglik = double(length=1),                                                               # --> double* loglik
                              weights = double(length=numstates),                                                      # --> double* weights
                              distr.type = as.integer(state.distributions),                                            # --> int* distr_type
                              size.initial = as.vector(size.initial),                                                  # --> double* initial_size
                              prob.initial = as.vector(prob.initial),                                                  # --> double* initial_prob
                              A.initial = as.vector(A.initial),                                                        # --> double* initial_A
                              proba.initial = as.vector(proba.initial),                                                # --> double* initial_proba
                              use.initial.params = as.logical(1),                                                      # --> bool* use_initial_params
                              num.threads = as.integer(num.threads),                                                   # --> int* num_threads
                              error = as.integer(0),                                                                   # --> int* error (error handling)
                              count.cutoff = as.integer(count.cutoff),                                                 # --> int* count.cutoff
                              algorithm = as.integer(algorithm),                                                       # --> int* algorithm
                              verbosity = as.integer(verbosity),                                                       # --> int* verbosity
                              PACKAGE = 'AneuFinder')
    hmm$eps             <- eps.try
    if(num.trials > 1){
      if(hmm$loglik.delta > hmm$eps){
        warning(paste0("ID = ",ID,": HMM did not converge in trial run ",i_try,"!\n"))
      }
      modellist[[as.character(i_try)]] <- hmm                                                                          # --> Store model in list
      init              <- 'random'
    }else if(num.trials == 1){
      if(hmm$loglik.delta > hmm$eps & istep == 1){
        warning(paste0("ID = ",ID,": HMM did not converge!\n"))
      }
    }                                                                                                                  # ==============================================================================
  }                                                                                                                    # RERUN MODEL WITH DIFFERENT INITIAL PARAMETERS
  if(num.trials > 1){                                                                                                  # ==============================================================================
    logliks             <- sapply(modellist,'[[','loglik')                                                             # --> Mathematically we should select the fit with highest loglikelihood. If we think the fit with the highest loglikelihood is incorrect, we should change the underlying model.
    df.weight           <- as.data.frame(lapply(modellist,'[[','weights'))                                             #     However, this is very complex and we choose to select a fit that we think is (more) correct, although it has not the highest support given our (imperfect) model.
    rownames(df.weight) <- state.labels
    models2use          <- (df.weight[most.frequent.state,] / apply(df.weight,2,max)) > 0.5                            # --> Select models where weight of most.frequent.state is at least half of that of actual most frequent state, then select model with highest loglik
    models2use[is.na(models2use)] <- FALSE                                                                             # --> Needed??????????????????????????
    if(any(models2use)){
      index2use         <- names(which.max(logliks[models2use]))
    }else{
      index2use         <- names(which.max(logliks))
    }
    hmm                 <- modellist[[index2use]]
    if(any(is.na(hmm$size) | is.nan(hmm$size) | is.infinite(hmm$size) |                                                # --> Check if size and prob parameter are correct
           is.na(hmm$prob) | is.nan(hmm$prob) | is.infinite(hmm$prob))){
      stop("...")
    }
    message(paste0("Rerunning trial ",index2use," with eps = ",eps))                                                   # --> Rerun the HMM with different epsilon and initial parameters from trial run
    hmm                 <- .C("C_univariate_hmm",
                              counts = as.integer(counts),                                                             # --> int* O
                              num.bins = as.integer(numbins),                                                          # --> int* T
                              num.states = as.integer(numstates),                                                      # --> int* N
                              state.labels = as.integer(state.labels),                                                 # --> int* state_labels
                              size = double(length=numstates),                                                         # --> double* size
                              prob = double(length=numstates),                                                         # --> double* prob
                              num.iterations = as.integer(max.iter),                                                   # --> int* maxiter
                              time.sec = as.integer(max.time),                                                         # --> double* maxtime
                              loglik.delta = as.double(eps),                                                           # --> double* eps
                              maxPosterior = double(length=numbins),                                                   # --> double* maxPosterior
                              states = integer(length=numbins),                                                        # --> int* states
                              A = double(length=numstates*numstates),                                                  # --> double* A
                              proba = double(length=numstates),                                                        # --> double* proba
                              loglik = double(length=1),                                                               # --> double* loglik
                              weights = double(length=numstates),                                                      # --> double* weights
                              distr.type = as.integer(state.distributions),                                            # --> int* distr_type
                              size.initial = as.vector(hmm$size),                                                      # --> double* initial_size
                              prob.initial = as.vector(hmm$prob),                                                      # --> double* initial_prob
                              A.initial = as.vector(hmm$A),                                                            # --> double* initial_A
                              proba.initial = as.vector(hmm$proba),                                                    # --> double* initial_proba
                              use.initial.params = as.logical(1),                                                      # --> bool* use_initial_params
                              num.threads = as.integer(num.threads),                                                   # --> int* num_threads
                              error = as.integer(0),                                                                   # --> int* error (error handling)
                              count.cutoff = as.integer(count.cutoff),                                                 # --> int* count.cutoff
                              algorithm = as.integer(algorithm),                                                       # --> int* algorithm
                              verbosity = as.integer(verbosity),                                                       # --> int* verbosity
                              PACKAGE = 'AneuFinder')                                                                  # ==============================================================================
  }                                                                                                                    # MAKE RETURN OBJECT
  result                <- list()                                                                                      # ==============================================================================
  result$ID             <- ID                                                                                          # --> ID
  result$bins           <- binned.data                                                                                 # --> Bin coordinates, counts and states
  result$bins$state     <- state.labels[hmm$states]
  result$bins$copy.number <- multiplicity[as.character(result$bins$state)]
  suppressMessages(                                                                                                    # --> Segmentation
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='copy.number',
                              columns2drop='width',columns2average=c('counts','mcounts','pcounts')),'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
  result$convergenceInfo <- list(eps=eps,loglik=hmm$loglik,loglik.delta=hmm$loglik.delta,                              # --> Convergence info
                              num.iterations=hmm$num.iterations,time.sec=hmm$time.sec,error=hmm$error)
  result$weights        <- hmm$weights                                                                                 # --> Weights
  names(result$weights) <- state.labels
  result$startProbs     <- hmm$proba                                                                                   # --> Probs
  names(result$startProbs) <- state.labels
  result$startProbs.initial <- hmm$proba.initial
  names(result$startProbs.initial) <- state.labels
  result$transitionProbs <- matrix(hmm$A,ncol=hmm$num.states)                                                          # --> Transition matrices
  rownames(result$transitionProbs) <- state.labels
  colnames(result$transitionProbs) <- state.labels
  result$transitionProbs.initial <- matrix(hmm$A.initial,ncol=hmm$num.states)
  rownames(result$transitionProbs.initial) <- state.labels
  colnames(result$transitionProbs.initial) <- state.labels
  distributions         <- data.frame()                                                                                # --> Distributions
  distributions.ini     <- data.frame()
  for(id in 1:length(hmm$distr.type)){
    distr               <- levels(state.distributions)[hmm$distr.type[id]]
    if(distr == 'dnbinom'){
      distributions     <- rbind(distributions,data.frame(type=distr,size=hmm$size[id],prob=hmm$prob[id],
                             mu=dnbinom.mean(hmm$size[id],hmm$prob[id]),
                             variance=dnbinom.variance(hmm$size[id],hmm$prob[id])))
      distributions.ini <- rbind(distributions.ini,data.frame(type=distr,size=hmm$size.initial[id],
                             prob=hmm$prob.initial[id],mu=dnbinom.mean(hmm$size.initial[id],hmm$prob.initial[id]),
                             variance=dnbinom.variance(hmm$size.initial[id],hmm$prob.initial[id])))
    }else if(distr == 'dgeom'){
      distributions     <- rbind(distributions,data.frame(type=distr,size=NA,prob=hmm$prob[id],
                             mu=dgeom.mean(hmm$prob[id]),variance=dgeom.variance(hmm$prob[id])))
      distributions.ini <- rbind(distributions.ini,data.frame(type=distr,size=NA,prob=hmm$prob.initial[id],
                             mu=dgeom.mean(hmm$prob.initial[id]),variance=dgeom.variance(hmm$prob.initial[id])))
    }else if (distr == 'delta'){
      distributions     <- rbind(distributions,data.frame(type=distr,size=NA,prob=NA,mu=0,variance=0))
      distributions.ini <- rbind(distributions.ini,data.frame(type=distr,size=NA,prob=NA,mu=0,variance=0))
    }else if (distr == 'dbinom'){
      distributions     <- rbind(distributions,data.frame(type=distr,size=hmm$size[id],prob=hmm$prob[id],
                             mu=dbinom.mean(hmm$size[id],hmm$prob[id]),
                             variance=dbinom.variance(hmm$size[id],hmm$prob[id])))
      distributions.ini <- rbind(distributions.ini,data.frame(type=distr,size=hmm$size.initial[id],
                             prob=hmm$prob.initial[id],mu=dbinom.mean(hmm$size.initial[id],hmm$prob.initial[id]),
                             variance=dbinom.variance(hmm$size.initial[id],hmm$prob.initial[id])))
    }
  }
  rownames(distributions) <- state.labels
  rownames(distributions.ini) <- state.labels
  result$distributions  <- distributions
  result$distributions.initial <- distributions.ini  
  result$warnings       <- warlist                                                                                     # --> Warnings and class
  class(result)         <- "aneuHMM"
  return(result)
}



RW_biHMM.findCNVs <- function(binned.data, ID=NULL, eps=0.01, init="standard", max.time=-1, max.iter=-1, num.trials=1, # ==============================================================================
                       eps.try=NULL, num.threads=1, count.cutoff.quantile=0.999,                                       # GET THE COUNTS AND CHECK THE DATA
                       states=c("zero-inflation",paste0(0:10,"-somy")), most.frequent.state="1-somy", algorithm='EM',  # ============================================================================== 
                       initial.params=NULL, verbosity=1){
  on.exit(.C("C_multivariate_cleanup", as.integer(num.comb.states), PACKAGE = 'AneuFinder'))                           # --> Define cleanup behaviour
  warlist               <- list()
  if(!is.null(initial.params)){
    init                <- 'initial.params'
  }
  if(is.null(eps.try)){
    eps.try             <- eps
  }
  counts                <- matrix(c(mcols(binned.data)[,'mcounts'],mcols(binned.data)[,'pcounts']),ncol=2,             # --> Get counts
                             dimnames=list(bin=1:length(binned.data), strand=c('minus','plus')))
  count.cutoff          <- ceiling(quantile(counts,count.cutoff.quantile))
  rows.higher           <- which(counts > count.cutoff)
  counts[rows.higher]   <- count.cutoff 
  if(length(rows.higher) > 0){
    message(paste0("Replaced read counts > ",count.cutoff," (quantile ",count.cutoff.quantile,") by ",count.cutoff,
	" in ",length(rows.higher)," bins. Set option 'count.cutoff.quantile=1' to disable this filtering.\n",
    "This filtering was done to enhance performance."))
  }
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data, otherwise HMM will blow up
    stop(paste0("ID = ",ID,": NAs found in reads."))
  }else if(all(counts == 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No HMM done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No HMM done."))
    result$warnings     <- warlist
    return(result)
  }                                                                                                                    # ==============================================================================
  if(init == 'initial.params'){                                                                                        # RUN UNIVARIATE HMM WHEN NO INITIAL PARAMETERS ARE GIVEN
    uni.transitionProbs <- initial.params$univariateParams$transitionProbs                                             # ==============================================================================
    uni.startProbs      <- initial.params$univariateParams$startProbs
    distributions       <- initial.params$distributions
    uni.weights         <- initial.params$univariateParams$weights
    uni.states          <- names(uni.weights)
    num.uni.states      <- length(uni.states)
    num.models          <- length(distributions)
    comb.states         <- factor(names(initial.params$startProbs),levels=names(initial.params$startProbs))
    num.comb.states     <- length(comb.states)
    states.list         <- list(minus=initial.params$bins$mstate, plus=initial.params$bins$pstate)
    comb.states.per.bin <- factor(do.call(paste,states.list),levels=levels(comb.states))
    A.initial           <- initial.params$transitionProbs
    proba.initial       <- initial.params$startProbs
    use.initial         <- TRUE
  }else if(init == 'standard'){
    binned.data.minus   <- binned.data                                                                                 # --> Stack the strands and run HMM.findCNVs
    strand(binned.data.minus) <- '-'
    binned.data.minus$counts <- binned.data.minus$mcounts
    binned.data.plus    <- binned.data
    strand(binned.data.plus) <- '+'
    binned.data.plus$counts <- binned.data.plus$pcounts
    binned.data.stacked <- c(binned.data.minus,binned.data.plus)
    message("Running univariate HMM")
    model.stacked       <- RW_HMM.findCNVs(binned.data.stacked,ID,eps=eps,init=init,max.time=max.time,                 # --> We here run an HMM on two sets of bins after eachother. A bit weird?????????
                             max.iter=max.iter,num.trials=num.trials,eps.try=eps.try,num.threads=num.threads,
                             count.cutoff.quantile=1,states=states,most.frequent.state=most.frequent.state)
    if(is.na(model.stacked$convergenceInfo$error)){                                                                    # --> Stop when univariate model gives errors.
      result$warnings   <- model.stacked$warnings                                                                      # --> What if not run? still an NA?? I think NULL
      return(result)
    }
    uni.transitionProbs <- model.stacked$transitionProbs                                                               # --> Extract probs, distributions
    uni.startProbs      <- model.stacked$startProbs
    distributions       <- model.stacked$distributions
    uni.weights         <- model.stacked$weights
    uni.states          <- names(uni.weights)
    num.uni.states      <- length(uni.states)
    num.models          <- 2
    levels.state        <- levels(model.stacked$bins$state)
    comb.states         <- paste(rep(levels.state,each=num.uni.states,times=1),
                             rep(levels.state,each=1,times=num.uni.states))
    comb.states         <- factor(comb.states, levels=comb.states)
    num.comb.states     <- length(comb.states)
    comb.states.per.bin <- paste(model.stacked$bins$state[as.character(strand(model.stacked$bins))=='-'],
                             model.stacked$bins$state[as.character(strand(model.stacked$bins))=='+'])
    comb.states.per.bin <- factor(comb.states.per.bin,levels=comb.states)
    A.initial           <- double(length=num.comb.states^2)
    proba.initial       <- double(length=num.comb.states)
    use.initial         <- FALSE                                                                                       # ==============================================================================
  }                                                                                                                    # CALCULATE DENSITIES THAT ARE NEEDED AS INPUT FOR MULTIVARIATE HMM
  z.per.count           <- array(NA,dim=c(count.cutoff+1,num.uni.states),                                              # ==============================================================================
                                 dimnames=list(counts=0:count.cutoff,state=uni.states))                                # --> Pre-compute z-values for each number of counts
  xcounts               <- 0:count.cutoff
  for(istate in 1:num.uni.states){
    if(distributions[istate,'type'] == 'dnbinom'){
      size              <- distributions[istate,'size']
      prob              <- distributions[istate,'prob']
      u                 <- stats::pnbinom(xcounts,size,prob)
    }else if(distributions[istate,'type'] == 'dbinom'){
      size              <- distributions[istate,'size']
      prob              <- distributions[istate,'prob']
      u                 <- stats::pbinom(xcounts,size,prob)
    }else if(distributions[istate,'type'] == 'delta'){
      u                 <- rep(1, length(xcounts))
    }else if(distributions[istate,'type'] == 'dgeom'){
      prob              <- distributions[istate,'prob']
      u                 <- stats::pgeom(xcounts, prob)
    }
    qnorm_u             <- stats::qnorm(u)
    mask                <- qnorm_u==Inf                                                                                # --> What about -Inf??????
    qnorm_u[mask]       <- stats::qnorm(1-1e-16)                                                                       # --> Why???
    z.per.count[,istate] <- qnorm_u
  }
  num.bins              <- length(binned.data)
  z.per.bin             <- array(NA,dim=c(num.bins,num.models,num.uni.states),                                         # --> Compute the z matrix
                             dimnames=list(bin=1:num.bins,strand=c("minus","plus"),state=uni.states))
  for(istrand in 1:num.models){
    for(istate in 1:num.uni.states){
      z.per.bin[,istrand,istate] <- z.per.count[counts[,istrand]+1,istate]
    }
  }
  correlationMatrix     <- array(0,dim=c(num.models,num.models,num.comb.states),                                       # --> Calculate correlation matrix
                             dimnames=list(strand=c("minus","plus"),strand=c("minus","plus"),comb.state=comb.states))
  correlationMatrixInverse <- correlationMatrix
  det.val               <- rep(0,num.comb.states)
  names(det.val)        <- comb.states
  for(comb.state in comb.states){
    state               <- strsplit(comb.state,' ')[[1]]
    mask                <- which(comb.states.per.bin == comb.state)
    z.temp              <- matrix(NA,ncol=num.models,nrow=length(mask))                                                # --> Subselect z
    for(istrand in 1:num.models){
      z.temp[,istrand]  <- z.per.bin[mask,istrand,state[istrand]]
    }
    temp <- tryCatch({
      if(nrow(z.temp) > 1){
        correlationMatrix[,,comb.state] <- cor(z.temp)
        det.val[comb.state] <- det(correlationMatrix[,,comb.state])
        correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
      }else{
        correlationMatrix[,,comb.state] <- diag(num.models)
        det.val[comb.state] <- det(correlationMatrix[,,comb.state])
        correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
      }
      0
    }, warning = function(war){
      1
    }, error = function(err){
      1
    })
    if(temp != 0){                                                                                                     # --> In case of an error or warning
      correlationMatrix[,,comb.state] <- diag(num.models)                                                              # --> Error or warning; ok.. but why did it occur?????
      det.val[comb.state] <- det(correlationMatrix[,,comb.state])
      correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
    }
  }
  densities             <- matrix(1,ncol=num.comb.states,nrow=num.bins,                                                # --> Calculate multivariate densities for each state
                                  dimnames=list(bin=1:num.bins,comb.state=comb.states))                                # --> At the end density has very small values!!!!! even e-300!!!!! Does this ever go wrong???
  for(comb.state in comb.states){
    istate              <- which(comb.state==comb.states)
    state               <- strsplit(comb.state,' ')[[1]]
    z.temp              <- matrix(NA,ncol=num.models,nrow=num.bins)
    product             <- 1
    for(istrand in 1:num.models){
      z.temp[,istrand]  <- z.per.bin[,istrand,state[istrand]]
      if(distributions[state[istrand],'type'] == 'dnbinom'){
        size            <- distributions[state[istrand],'size']
        prob            <- distributions[state[istrand],'prob']
        product         <- product * stats::dnbinom(counts[,istrand],size,prob)
      }else if(distributions[state[istrand],'type'] == 'dgeom'){
        prob            <- distributions[state[istrand],'prob']
        product         <- product * stats::dgeom(counts[,istrand],prob)
      }else if(distributions[state[istrand],'type'] == 'delta'){
        product         <- product * ifelse(counts[,istrand]==0, 1, 0)
      }
    }
    exponent            <- -0.5 * apply((z.temp %*% (correlationMatrixInverse[,,istate]-diag(num.models)))*z.temp,
                             1,sum)
    exponent[is.nan(exponent)] <- 0
    densities[,istate]  <- product * det.val[istate]^(-0.5) * exp(exponent)
  }
  densities[densities > 1] <- 1                                                                                        # --> Check if densities are > 1 or < 0
  densities[densities < 0] <- 0                                                                                        # --> Is it ok if this happens???????
  check                 <- which(apply(densities, 1, sum) == 0)                                                        # --> Check if densities are 0 everywhere in some bins
  if(length(check) > 0){
    if(check[1] == 1){
      densities[1,]     <- rep(1e-10, ncol(densities))                                                                 # --> Why 10**-10??????
      check             <- check[-1]
    }
    for(icheck in check){
      densities[icheck,] <- densities[icheck-1,]                                                                       # --> If row sum is zero then take values of previous row. Why?????????
    }
  }                                                                                                                    # ==============================================================================
  algorithm             <- factor(algorithm, levels=c('baumWelch','viterbi','EM'))                                     # RUN THE MULTIVARIATE HMM
  hmm                   <- .C("C_multivariate_hmm",                                                                    # ==============================================================================
                              densities = as.double(densities),                                                        # --> double* D
                              num.bins = as.integer(num.bins),                                                         # --> int* T
                              num.comb.states = as.integer(num.comb.states),                                           # --> int* N
                              num.strands = as.integer(num.models),                                                    # --> int* Nmod
                              comb.states = as.integer(comb.states),                                                   # --> int* comb_states
                              num.iterations = as.integer(max.iter),                                                   # --> int* maxiter
                              time.sec = as.integer(max.time),                                                         # --> double* maxtime
                              loglik.delta = as.double(eps),                                                           # --> double* eps
                              maxPosterior = double(length=num.bins),                                                  # --> double* maxPosterior
                              states = integer(length=num.bins),                                                       # --> int* states
                              A = double(length=num.comb.states*num.comb.states),                                      # --> double* A
                              proba = double(length=num.comb.states),                                                  # --> double* proba
                              loglik = double(length=1),                                                               # --> double* loglik
                              A.initial = as.vector(A.initial),                                                        # --> double* initial_A
                              proba.initial = as.vector(proba.initial),                                                # --> double* initial_proba
                              use.initial.params = as.logical(use.initial),                                            # --> bool* use_initial_params
                              num.threads = as.integer(num.threads),                                                   # --> int* num_threads
                              error = as.integer(0),                                                                   # --> error handling
                              algorithm = as.integer(algorithm),                                                       # --> int* algorithm
                              verbosity = as.integer(verbosity),                                                       # --> int* verbosity
                              PACKAGE = 'AneuFinder')
  if(hmm$loglik.delta > eps){                                                                                          # --> Check convergence
    warning(paste0("ID = ",ID,": HMM did not converge!\n"))
  }
  if(hmm$error == 1){
    warlist[[length(warlist)+1]] <- 
      warning(paste0("ID = ",ID,": A NaN occurred during the Baum-Welch! Parameter estimation terminated prematurely.",
      "\n Check your library! The following factors are known to cause this error:\n",
      "1) Your read counts contain very high numbers. Try again with a lower value for 'count.cutoff.quantile'.\n",
      "2) Your library contains too few reads in each bin.\n",
      "3) Your library contains reads for a different genome than it was aligned to."))
    result$warnings     <- warlist
    return(result)
  }else if(hmm$error == 2){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": An error occurred during the Baum-Welch!\n",
      "Parameter estimation terminated prematurely. Check your library."))
    result$warnings     <- warlist
    return(result)
  }                                                                                                                    # ==============================================================================
  result                <- list()                                                                                      # MAKE RETURN OBJECT
  result$ID             <- ID                                                                                          # ==============================================================================
  result$bins           <- binned.data                                                                                 # --> Bin coordinates and counts
  matrix.states         <- do.call("rbind",strsplit(as.character(comb.states[hmm$states]),split=' '))                  # --> Get states as factors in data.frame
  mstate.num            <- sub("zero-inflation","0-somy",matrix.states[,1])
  mstate.num            <- as.numeric(sub("-somy", "", mstate.num))
  pstate.num            <- sub("zero-inflation","0-somy",matrix.states[,2])
  pstate.num            <- as.numeric(sub("-somy", "", pstate.num))
  bs.copy.num           <- mstate.num + pstate.num
  bs.state              <- paste0(bs.copy.num,"-somy")
  result$bins$state     <- factor(bs.state, levels=unique(c(states,sort(unique(bs.state)))))
  result$bins$mstate    <- factor(matrix.states[,1], levels=uni.states)
  result$bins$pstate    <- factor(matrix.states[,2], levels=uni.states)
  result$bins$copy.number <- bs.copy.num
  result$bins$mcopy.number <- mstate.num
  result$bins$pcopy.number <- pstate.num      	
  result$bins$combi     <- paste(result$bins$mcopy.number, result$bins$pcopy.number)                                   # --> Segmentation
  suppressMessages(
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='combi',
                            columns2drop='width',columns2average=c('counts','mcounts','pcounts')),'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]
  bin.num               <- NULL                                                                                        # --> Determine the number of bins for each segment
  for(chr in unique(as.character(seqnames(result$bins)))){
    bin.num             <- c(bin.num,rle(result$bins$combi[which(as.character(seqnames(result$bins)) == chr)])$lengths)
  }
  result$segments$num.bins <- bin.num
  mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts","mean.pcounts","state",   # --> Re-order the columns
                                                     "mstate","pstate","copy.number","mcopy.number","pcopy.number")]
  result$bins$combi     <- NULL
  result$convergenceInfo <- list(eps=eps,loglik=hmm$loglik,loglik.delta=hmm$loglik.delta,                              # --> Convergence info
                                 num.iterations=hmm$num.iterations,time.sec=hmm$time.sec)
  result$weights        <- table(result$bins$state)/length(result$bins)                                                # --> Weights
  result$startProbs     <- hmm$proba                                                                                   # --> Initial probs
  names(result$startProbs) <- comb.states
  result$startProbs.initial <- hmm$proba.initial
  names(result$startProbs.initial) <- comb.states
  result$transitionProbs <- matrix(hmm$A, ncol=num.comb.states)                                                        # --> Transition matrices
  colnames(result$transitionProbs) <- comb.states
  rownames(result$transitionProbs) <- comb.states
  result$transitionProbs.initial <- matrix(hmm$A.initial, ncol=num.comb.states)
  colnames(result$transitionProbs.initial) <- comb.states
  rownames(result$transitionProbs.initial) <- comb.states
  result$distributions  <- distributions                                                                               # --> Distributions
  result$univariateParams <- list(weights=uni.weights,startProbs=uni.startProbs,                                       # --> Univariate infos
                                  transitionProbs=uni.transitionProbs,distributions=distributions)
  result$warnings       <- warlist
  class(result)         <- "aneuBiHMM"
  return(result)
}



RW_DNAcopy.findCNVs <- function(binned.data, ID=NULL, CNgrid.start=1.5, strand='*'){                                   # ==============================================================================
  warlist               <- list()                                                                                      # GET THE COUNTS AND CHECK THE DATA
  if(strand=='+'){                                                                                                     # ==============================================================================
    select              <- 'pcounts'                                                                                   # --> Get counts
  }else if(strand=='-'){
    select              <- 'mcounts'
  }else if(strand=='*'){
    select              <- 'counts'
  }
  counts                <- mcols(binned.data)[,select]
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data
    stop(paste0("ID = ",ID,": NAs found in reads."))                                                                   # --> Also important for DNAcopy?????
  }else if(all(counts == 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No DNAcopy done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No DNAcopy done."))
    result$warnings     <- warlist                                                                                     # ==============================================================================
    return(result)                                                                                                     # RUN DNACOPY
  }                                                                                                                    # ==============================================================================
  set.seed(0)                                                                                                          # --> Fix seed to get reproducible results
  counts.normal         <- (counts+1) / mean(counts+1)
  logcounts             <- log2(counts.normal)
  CNA.object            <- DNAcopy::CNA(genomdat=logcounts,maploc=as.numeric(start(binned.data)),data.type='logratio',
                             chrom=paste0(as.character(seqnames(binned.data)),"_",as.character(strand(binned.data))))  # --> We need strand information for the bivariate variant of DNA copy (stacked strands).
  CNA.smoothed          <- DNAcopy::smooth.CNA(CNA.object)
  CNA.segs              <- DNAcopy::segment(CNA.smoothed,verbose=0,min.width=5)
  CNA.segs              <- CNA.segs$output
  CNA.segs$chrom        <- as.character(CNA.segs$chrom)
  CNA.segs$strand       <- substr(CNA.segs$chrom,nchar(CNA.segs$chrom),nchar(CNA.segs$chrom))
  CNA.segs$chrom        <- substr(CNA.segs$chrom,1,nchar(CNA.segs$chrom)-2)
  CNA.segs$loc.end      <- c(CNA.segs$loc.start[2:nrow(CNA.segs)] - 1, NA)                                             # --> End position is start next - 1
  last.ind.chr.strand   <- cumsum(rle(paste0(CNA.segs$chrom,"_",CNA.segs$strand))$lengths)
  CNA.segs$loc.end[last.ind.chr.strand] <- seqlengths(binned.data)[CNA.segs$chrom[last.ind.chr.strand]]                # --> End position last segment of each chromosome is the length of the chromosome
  segs.gr               <- GRanges(seqnames=CNA.segs$chrom,
                             ranges=IRanges(start=CNA.segs$loc.start,end=CNA.segs$loc.end),strand=CNA.segs$strand)
  segs.gr$mean.count    <- (2^CNA.segs$seg.mean) * mean(counts+1) - 1
  segs.gr$mean.count[segs.gr$mean.count < 0] <- 0
  seg.ind               <- findOverlaps(binned.data,segs.gr,select='first')                                            # --> Modify bins to contain median count
  counts.normal         <- counts / mean(counts[which(counts > 0)])
  segs.gr$median.count  <- sapply(split(counts.normal,seg.ind),function(x){                                            # --> Bit ugly approach!!! And it is definitely not the median!!!
    qus                 <- quantile(x, c(0.01, 0.99))                                                                  # --> Why not the median?? I guess you want to get rid of outliers when possible??
    y                   <- x[x >= qus[1] & x <= qus[2]]
    if(sum(y) == 0 | length(y) == 0){
      y                 <- x
    }
    mu                  <- mean(y)
    return(mu)
  })
  counts.median         <- segs.gr$median.count[seg.ind]
  CNgrid                <- seq(CNgrid.start,6,by=0.01)                                                                 # --> Determine copy number
  outerRaw              <- counts.median %o% CNgrid                                                                    # --> Multiplication of two vectors. First with first, first with second, etc. 
  outerDiff             <- (outerRaw - round(outerRaw)) ^ 2
  sumOfSquares          <- colSums(outerDiff,na.rm=FALSE,dims=1)                                                       # --> For each muliplication factor the sum of squared differences
  CN                    <- CNgrid[order(sumOfSquares)[1]]                                                              # --> Pick the factor that gave the smallest sum of squared differences.  ??????? Is sum of squared differences enough????? For edivisve we also use multiplication factor and the number of segments....!!!!
  CN.states             <- round(counts.median * CN)                                                                   # --> Determine states
  names(CN.states)      <- NULL                                                                                        # ==============================================================================
  result                <- list()                                                                                      # MAKE RETURN OBJECT
  result$ID             <- ID                                                                                          # ==============================================================================
  result$bins           <- binned.data                                                                                 # --> Bin coordinates and counts
  result$bins$state     <- factor(paste0(CN.states, '-somy'),levels=paste0(sort(unique(CN.states)),'-somy'))
  result$bins$copy.number <- CN.states
  suppressMessages(                                                                                                    # --> Segmentation
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='state',
                             columns2drop='width',columns2average=c('counts','mcounts','pcounts')), 'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
  chr.strand            <- paste0(as.character(seqnames(result$bins)),"_",as.character(strand(result$bins)))           # --> Determine the number of bins for each segment
  bin.num               <- NULL
  for(chr.st in unique(chr.strand)){
    bin.num             <- c(bin.num,rle(as.character(result$bins$state[which(chr.strand == chr.st)]))$lengths)
  }
  result$segments$num.bins <- bin.num
  mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts",                          # --> Re-order the columns
                                                     "mean.pcounts","state","copy.number")]
  result$weights        <- table(result$bins$state) / length(result$bins)                                              # --> Weights
  result$distributions  <- assign.distributions(counts=result$bins$counts,states=result$bins$state)                    # --> Distributions
  result$warnings       <- warlist
  class(result)         <- "aneuHMM"                                                                                   # --> Ugly!!!!!
  return(result)  
}



RW_biDNAcopy.findCNVs <- function(binned.data,ID=NULL,CNgrid.start=0.5){                                               # ==============================================================================
  warlist               <- list()                                                                                      # GET AND CHECK THE COUNTS
  counts                <- matrix(c(mcols(binned.data)[,'mcounts'],mcols(binned.data)[,'pcounts']),ncol=2,             # ==============================================================================
                             dimnames=list(bin=1:length(binned.data),strand=c('minus','plus')))                        # --> Get counts
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data
    stop(paste0("ID = ",ID,": NAs found in reads."))                                                                   # --> Also important for DNAcopy?????
  }else if(all(counts == 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No DNAcopy done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No DNAcopy done."))
    result$warnings     <- warlist
    return(result)
  }                                                                                                                    # ==============================================================================
  binned.data.minus     <- binned.data                                                                                 # STACK THE STRANDS AND RUN DNACOPY
  strand(binned.data.minus) <- '-'                                                                                     # ==============================================================================
  binned.data.minus$counts <- binned.data.minus$mcounts
  binned.data.plus      <- binned.data
  strand(binned.data.plus) <- '+'
  binned.data.plus$counts <- binned.data.plus$pcounts
  binned.data.stacked   <- c(binned.data.minus,binned.data.plus)
  message("Running DNAcopy")
  model.stacked         <- RW_DNAcopy.findCNVs(binned.data.stacked,ID,CNgrid.start=CNgrid.start)                       # ==============================================================================
  result                <- list()                                                                                      # MAKE RETURN OBJECT
  result$ID             <- ID                                                                                          # ==============================================================================
  result$bins           <- binned.data
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
  result$bins$combi     <- paste(result$bins$mcopy.number,result$bins$pcopy.number)                                    # --> Segmentation
  suppressMessages(
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='combi',
                            columns2drop='width',columns2average=c('counts','mcounts','pcounts')),'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]
  bin.num               <- NULL                                                                                        # --> Determine the number of bins for each segment
  for(chr in unique(as.character(seqnames(result$bins)))){
    bin.num             <- c(bin.num,rle(result$bins$combi[which(as.character(seqnames(result$bins)) == chr)])$lengths)
  }
  result$segments$num.bins <- bin.num
  mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts","mean.pcounts","state",   # --> Re-order the columns
                                                     "mstate","pstate","copy.number","mcopy.number","pcopy.number")]
  result$bins$combi     <- NULL
  result$weights        <- table(result$bins$state) / length(result$bins)                                              # --> Weights
  result$distributions  <- list(minus = model.stacked$distributions, plus = model.stacked$distributions)               # --> Distributions for strands separate
  result$distributions$both <- assign.distributions(counts=result$bins$counts,states=result$bins$state)                # --> Distributions for strands combined
  result$univariateParams <- list(weights=model.stacked$weights)
  result$warnings       <- model.stacked$warnings
  class(result)         <- "aneuBiHMM"                                                                                 # --> Ugly!!!!!
  return(result)
}



assign.distributions <- function(counts,states){
  distributions         <- list()
  bins.splt             <- split(counts,states)                                                                        # --> Determine distribution type based on mu and variance.
  for(i1 in 1:length(bins.splt)){                                                                                      # --> Is this correct?????? Why do we need this? Originally the results were not based on these distributions, right?
    if(length(bins.splt[[i1]]) > 0){
      qus               <- quantile(bins.splt[[i1]], c(0.01, 0.99))
      qcounts           <- bins.splt[[i1]]
      qcounts           <- qcounts[qcounts >= qus[1] & qcounts <= qus[2]]
      if(sum(qcounts) == 0 | length(qcounts)==0){
        qcounts         <- bins.splt[[i1]]
      }
      mu                <- mean(qcounts)
      variance          <- var(qcounts)
      if(is.na(variance)){
        variance        <- mu + 1                                                                                      # --> Somewhat arbitrary. A bit strange?????
      }
      if(names(bins.splt)[i1] == '0-somy'){
        distr           <- 'dgeom'
        size            <- NA
        prob            <- dgeom.prob(mu)
      }else{
        if(is.na(variance) | is.na(mu)){
          distr         <- 'dnbinom'
          size          <- NA
          prob          <- NA
        }else{
          if(variance < mu){
            distr       <- 'dbinom'
            size        <- dbinom.size(mu,variance)
            prob        <- dbinom.prob(mu,variance)
          }else if(variance > mu){
            distr       <- 'dnbinom'
            size        <- dnbinom.size(mu,variance)
            prob        <- dnbinom.prob(mu,variance)
          }else{
            distr       <- 'dpois'
            size        <- NA
            prob        <- mu
          }
        }
      }
      distributions[[i1]] <- data.frame(type=distr,size=size,prob=prob,mu=mu,variance=variance)
    }
  }
  distributions         <- do.call(rbind,distributions)
  rownames(distributions) <- names(which(lapply(bins.splt,length) > 0))
  return(distributions)
}



RW_edivisive.findCNVs <- function(binned.data,ID=NULL,CNgrid.start=1.5,strand='*',R=10,sig.lvl=0.1){                   # ==============================================================================
  warlist               <- list()                                                                                      # GET AND CHECK THE COUNTS
  if(strand=='+'){                                                                                                     # ==============================================================================
    select              <- 'pcounts'                                                                                   # --> Get counts
  }else if(strand=='-'){
    select              <- 'mcounts'
  }else if(strand=='*'){
    select              <- 'counts'
  }
  counts                <- mcols(binned.data)[,select]
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data
    stop(paste0("ID = ",ID,": NAs found in reads."))                                                                   # --> Also important for Edivisive?????
  }else if(all(counts == 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No Edivisive done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No Edivisive done."))
    result$warnings     <- warlist                                                                                     # ==============================================================================
    return(result)                                                                                                     # 
  }                                                                                                                    # ==============================================================================
  set.seed(0)                                                                                                          # --> Fix seed to get reproducible results
  binned.data$cluster   <- NA
  cl                    <- 0
  for(chrom in unique(as.character(seqnames(binned.data)))){
    chr.rows            <- which(as.character(seqnames(binned.data)) == chrom)
    counts.chrom        <- counts[chr.rows]
    dim(counts.chrom)   <- c(length(counts.chrom),1)
    cp                  <- ecp::e.divisive(counts.chrom,min.size=2,R=R,sig.lvl=sig.lvl)
    binned.data$cluster[chr.rows] <- cp$cluster + cl
    cl                  <- cl + length(cp$p.values)
  }
  counts.normal         <- counts / mean(counts[which(counts > 0)])                                                    # --> Modify bins to contain mean count
  cnmean                <- sapply(split(counts.normal,binned.data$cluster),function(x){
    qus                 <- quantile(x, c(0.01, 0.99))
    y                   <- x[x >= qus[1] & x <= qus[2]]
    if(sum(y) == 0 | length(y) == 0){
      y                 <- x
    }
    mu                  <- mean(y)
    return(mu)
  })
  counts.normal.mean    <- cnmean[as.character(binned.data$cluster)]
  CNgrid                <- seq(CNgrid.start,6,by=0.01)                                                                 # --> Determine copy number
  outerRaw              <- counts.normal.mean %o% CNgrid
  outerDiff             <- (outerRaw - round(outerRaw)) ^ 2
  sumOfSquares          <- colSums(outerDiff,na.rm=FALSE,dims=1)
  CN                    <- CNgrid[order(sumOfSquares)][1]
  CN.states             <- round(counts.normal.mean * CN)
  names(CN.states)      <- NULL                                                                                        # ==============================================================================
  result                <- list()                                                                                      # MAKE RETURN OBJECT
  result$ID             <- ID                                                                                          # ==============================================================================
  result$bins           <- binned.data
  result$bins$state     <- factor(paste0(CN.states,'-somy'), levels=paste0(sort(unique(CN.states)),'-somy'))
  result$bins$copy.number <- CN.states
  suppressMessages(
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='state',
                             columns2drop='width',columns2average=c('counts','mcounts','pcounts')),'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
  bin.num               <- NULL                                                                                        # --> Determine the number of bins for each segment
  for(chr in unique(as.character(seqnames(result$bins)))){
    bin.num             <- c(bin.num,rle(as.character(result$bins$state[
                             which(as.character(seqnames(result$bins)) == chr)]))$lengths)
  }
  result$segments$num.bins <- bin.num
  mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts",                          # --> Re-order the columns
                                                     "mean.pcounts","state","copy.number")]
  result$bins$cluster   <- NULL
  result$weights        <- table(result$bins$state) / length(result$bins)                                              # --> Weights
  result$distributions  <- assign.distributions(counts=result$bins$counts,states=result$bins$state)                    # --> Distributions
  result$warnings       <- warlist
  class(result)         <- "aneuHMM"                                                                                   # --> Ugly!!!!!
  return(result)
}



RW_bi.edivisive.findCNVs <- function(binned.data,ID=NULL,CNgrid.start=0.5,R=10,sig.lvl=0.1){                           # ==============================================================================
  warlist               <- list()                                                                                      # GET AND CHECK THE COUNTS
  counts                <- as.matrix(mcols(binned.data)[,c('mcounts','pcounts')])                                      # ==============================================================================
  if(any(is.na(counts))){                                                                                              # --> Check if there are counts in the data  --> Check with other bi functions!!!
    stop(paste0("ID = ",ID,": NAs found in reads."))                                                                   # --> Also important for Edivisive?????
  }else if(all(counts == 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": All counts in data are zero. No Edivisive done."))
    result$warnings     <- warlist
    return(result)
  }else if (any(counts < 0)){
    warlist[[length(warlist)+1]] <- warning(paste0("ID = ",ID,": Some counts in data are negative. No Edivisive done."))
    result$warnings     <- warlist
    return(result)
  }                                                                                                                    # ==============================================================================
  set.seed(0)                                                                                                          # --> Fix seed to get reproducible results
  binned.data$cluster   <- NA
  cl                    <- 0
  for(chrom in unique(as.character(seqnames(binned.data)))){
    chr.rows            <- which(as.character(seqnames(binned.data)) == chrom)
    counts.chrom        <- counts[chr.rows,]
    cp                  <- ecp::e.divisive(counts.chrom,min.size=2,R=R,sig.lvl=sig.lvl)
    binned.data$cluster[chr.rows] <- cp$cluster + cl
    cl                  <- cl + length(cp$p.values)
  }
  counts.normal         <- counts / mean(counts[which(counts > 0)])                                                    # --> Modify bins to contain median count
  cnmean.m              <- numeric()
  cnmean.p              <- numeric()
  for(i1 in 1:max(binned.data$cluster)){
    x                   <- counts.normal[binned.data$cluster==i1,,drop=FALSE]
    qus                 <- quantile(x, c(0.01, 0.99))
    within.quantile     <- apply(x,2,function(z){z >= qus[1] & z <= qus[2]})                                           # --> The counts of both strands need to be within the quantile range
    dim(within.quantile) <- dim(x)                                                                                     # --> Keep the dimensions of thw matrix
    dimnames(within.quantile) <- dimnames(x)
    within.quantile     <- within.quantile[,'mcounts'] & within.quantile[,'pcounts']
    y                   <- x[within.quantile,,drop=FALSE]
    if(sum(y) == 0 | length(y)==0){
      y                 <- x
    }
    mu                  <- colMeans(y)
    cnmean.m[as.character(i1)] <- mu['mcounts']
    cnmean.p[as.character(i1)] <- mu['pcounts']
  }
  counts.normal.mean.m  <- cnmean.m[as.character(binned.data$cluster)]
  counts.normal.mean.p  <- cnmean.p[as.character(binned.data$cluster)]
  counts.normal.mean.stacked <- c(counts.normal.mean.m,counts.normal.mean.p)
  CNgrid                <- seq(CNgrid.start,6,by=0.01)                                                                 # --> Determine copy number
  outerRaw              <- counts.normal.mean.stacked %o% CNgrid
  outerDiff             <- (outerRaw - round(outerRaw)) ^ 2
  sumOfSquares          <- colSums(outerDiff,na.rm=FALSE,dims=1)
  CN                    <- CNgrid[order(sumOfSquares)][1]
  CN.states             <- round(counts.normal.mean.stacked * CN)
  names(CN.states)      <- NULL                                                                                        # ==============================================================================
  result                <- list()                                                                                      # MAKE RETURN OBJECT
  result$ID             <- ID                                                                                          # ==============================================================================
  result$bins           <- binned.data
  result$bins$mcopy.number <- CN.states[1:(length(CN.states)/2)]
  result$bins$pcopy.number <- CN.states[((length(CN.states)/2)+1):length(CN.states)]
  result$bins$copy.number  <- result$bins$mcopy.number + result$bins$pcopy.number
  lev                   <- paste0(0:max(result$bins$copy.number),'-somy')
  result$bins$state     <- factor(paste0(result$bins$copy.number,'-somy'),levels=lev)
  result$bins$mstate    <- factor(paste0(result$bins$mcopy.number,'-somy'),levels=lev)
  result$bins$pstate    <- factor(paste0(result$bins$pcopy.number,'-somy'),levels=lev)
  result$bins$combi     <- paste(result$bins$mcopy.number,result$bins$pcopy.number)                                    # --> Segmentation
  suppressMessages(
    result$segments     <- as(collapseBins(as.data.frame(result$bins),column2collapseBy='combi',
                             columns2drop='width',columns2average=c('counts','mcounts','pcounts')),'GRanges')
  )
  seqlevels(result$segments) <- seqlevels(result$bins)                                                                 # --> Correct order from as()
  seqlengths(result$segments) <- seqlengths(result$bins)[names(seqlengths(result$segments))]
  bin.num               <- NULL                                                                                        # --> Determine the number of bins for each segment
  for(chr in unique(as.character(seqnames(result$bins)))){
    bin.num             <- c(bin.num,rle(result$bins$combi[which(as.character(seqnames(result$bins)) == chr)])$lengths)
  }
  result$segments$num.bins <- bin.num
  mcols(result$segments) <- mcols(result$segments)[c("num.bins","mean.counts","mean.mcounts","mean.pcounts","state",   # --> Re-order the columns of segments
                                                     "mstate","pstate","copy.number","mcopy.number","pcopy.number")]
  mcols(result$bins)    <- mcols(result$bins)[c("counts","mcounts","pcounts","GC","state","mstate",                    # --> Re-order the columns of bins
                                                "pstate","copy.number","mcopy.number","pcopy.number")]
  result$weights        <- table(result$bins$state) / length(result$bins)                                              # --> Weights
  dist.strand           <- assign.distributions(counts=c(result$bins$mcounts,result$bins$pcounts),                     # --> Distributions for strand separately
                             states=c(as.character(result$bins$mstate),as.character(result$bins$pstate)))
  result$distributions  <- list(minus = dist.strand, plus = dist.strand)
  result$distributions$both <- assign.distributions(counts=result$bins$counts,states=result$bins$state)                # --> Distributions for strands combined
  result$univariateParams <- list(weights=table(c(as.character(result$bins$mstate),
                               as.character(result$bins$pstate))) / (length(result$bins)*2))
  result$warnings       <- warlist
  class(result)         <- "aneuBiHMM"                                                                                 # --> Ugly!!!!!
  return(result)
}





#binned.data <- get(load("D:\\Rene\\Projects\\AneuFinder\\Data\\Data_test_new_code\\Bam_out\\binned-GC\\MB180503_I_AD_001.bam.RData"))
#binned.data <- binned.data[[2]]

#ID                    <- "ABC"
#method                <- "HMM"
#strand                <- '*'
#R                     <- 10
#sig.lvl               <- 0.1
#eps                   <- 0.01
#init                  <- "standard"
#max.time              <- 60
#max.iter              <- -1
#num.trials            <- 15
#eps.try               <- max(10*eps, 1)
#num.threads           <- 1
#count.cutoff.quantile <- 0.999
#states                <- c("zero-inflation",paste0(0:10,"-somy"))
#most.frequent.state   <- "2-somy"
#algorithm             <- "EM"
#initial.params        <- NULL
#verbosity             <- 1

#CNgrid.start=0.5

#library(AneuFinder)



#plot(sumOfSquares,pch=20,cex=5,col=rgb(0,0.5,0,0.1))




