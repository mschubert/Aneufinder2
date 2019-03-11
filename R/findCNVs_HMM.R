#' Find copy number variations (univariate)
#'
#' \code{HMM.findCNVs} classifies the binned read counts into several states
#' which represent copy-number-variation.
#'
#' @param binned.data A \code{\link{GRanges-class}} object with binned read
#'   counts. Alternatively a \code{\link{GRangesList}} object with offsetted
#'   read counts.
#' @param ID An identifier that will be used to identify this sample in various
#'   downstream functions. Could be the file name of the \code{binned.data} for
#'   example.
#' @param eps method-HMM: Convergence threshold for the Baum-Welch algorithm.
#' @param init method-HMM: One of the following initialization procedures:
#'  \describe{
#'    \item{\code{standard}}{The negative binomial of state '2-somy' will be
#'      initialized with \code{mean=mean(counts)}, \code{var=var(counts)}. This
#'      procedure usually gives good convergence.}
#'    \item{\code{random}}{Mean and variance of the negative binomial of state
#'     '2-somy' will be initialized with random values (in certain boundaries,
#'     see source code). Try this if the \code{standard} procedure fails to
#'     produce a good fit.}
#'  }
#' @param max.time method-HMM: The maximum running time in seconds for the
#'   Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate
#'   after the current iteration finishes. Set \code{max.time = -1} for no limit.
#' @param max.iter method-HMM: The maximum number of iterations for the
#'   Baum-Welch algorithm. Set \code{max.iter = -1} for no limit.
#' @param num.trials method-HMM: The number of trials to find a fit where state
#'   \code{most.frequent.state} is most frequent. Each time, the HMM is seeded
#'   with different random initial values.
#' @param eps.try method-HMM: If code num.trials is set to greater than 1,
#'   \code{eps.try} is used for the trial runs. If unset, \code{eps} is used.
#' @param num.threads method-HMM: Number of threads to use. Setting this to >1
#'   may give increased performance.
#' @param count.cutoff.quantile method-HMM: A quantile between 0 and 1. Should
#'   be near 1. Read counts above this quantile will be set to the read count
#'   specified by this quantile. Filtering very high read counts increases the
#'   performance of the Baum-Welch fitting procedure. However, if your data
#'   contains very few peaks they might be filtered out. Set
#'   \code{count.cutoff.quantile=1} in this case.
#' @param strand Find copy-numbers only for the specified strand. One of
#'   \code{c('+', '-', '*')}.
#' @param states method-HMM: A subset or all of
#'   \code{c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy",...)}.
#'   This vector defines the states that are used in the Hidden Markov Model.
#'   The order of the entries must not be changed.
#' @param most.frequent.state method-HMM: One of the states that were given in
#'   \code{states}. The specified state is assumed to be the most frequent one.
#'   This can help the fitting procedure to converge into the correct fit.
#' @param algorithm method-HMM: One of \code{c('baumWelch','EM')}. The
#'   expectation maximization (\code{'EM'}) will find the most likely states and
#'   fit the best parameters to the data, the \code{'baumWelch'} will find the
#'   most likely states using the initial parameters.
#' @param initial.params method-HMM: A \code{\link{aneuHMM}} object or file
#'   containing such an object from which initial starting parameters will be
#'   extracted.
#' @param verbosity method-HMM: Integer specifying the verbosity of printed
#'   messages.
#' @return An \code{\link{aneuHMM}} object.
#' @importFrom stats runif
HMM.findCNVs <- function(binned.data, ID=NULL, eps=0.01, init="standard",
                         max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL,
                         num.threads=1, count.cutoff.quantile=0.999,
                         strand='*',
                         states=c("zero-inflation",paste0(0:10,"-somy")),
                         most.frequent.state="2-somy", algorithm="EM",
                         initial.params=NULL, verbosity=1) {

    on.exit(.C("C_univariate_cleanup", PACKAGE = 'AneuFinder'))

    warlist <- list()
    if (!is.null(initial.params))
        init <- 'initial.params'
    if (num.trials == 1)
        eps.try <- eps

    if (strand=='+')
        select <- 'pcounts'
    else if (strand=='-')
        select <- 'mcounts'
    else if (strand=='*')
        select <- 'counts'

    # --> Filter high counts out, makes HMM faster
    counts <- mcols(binned.data)[,select]
    count.cutoff <- ceiling(quantile(counts, count.cutoff.quantile))
    rows.higher           <- which(counts > count.cutoff)
    counts[rows.higher]   <- count.cutoff

    if (length(rows.higher) > 0)
        message("Replaced read counts > ",count.cutoff," (quantile ", count.cutoff.quantile,
                ") by ",count.cutoff, " in ", length(rows.higher),
                " bins. Set option 'count.cutoff.quantile=1' to disable this filtering.",
                " This filtering was done to enhance performance.")

    # --> Check if there are counts in the data, otherwise HMM will blow up
    if (any(is.na(counts))) {
        stop(paste0("ID = ",ID,": NAs found in reads."))
    } else if(!any(counts != 0)) {
        warlist[[length(warlist)+1]] <- warning("ID = ",ID,": All counts in data are zero. No HMM done.")
        result$warnings <- warlist
        return(result)
    } else if (any(counts < 0)) {
        warlist[[length(warlist)+1]] <- warning("ID = ",ID,": Some counts in data are negative. No HMM done.")
        result$warnings <- warlist
        return(result)
    }

    numbins <- length(binned.data)
    numstates <- length(states)

    inistates <- initializeStates(states)
    state.labels <- inistates$states
    state.distributions <- inistates$distributions
    multiplicity <- inistates$multiplicity
    dependent.states.mask <- (state.labels != 'zero-inflation') & (state.labels != '0-somy')

    # --> This factor is later on changed into a number.
    algorithm <- factor(algorithm,levels=c('baumWelch','viterbi','EM'))

    # RUN MODEL WITH EACH TIME DIFFERENT INITIAL PARAMETERS
    modellist <- list()
    for(i_try in 1:num.trials) {
        if (verbosity >= 1)
            message(paste0("Trial ",i_try," / ",num.trials))

        if(init == 'initial.params') {
            A.initial <- initial.params$transitionProbs
            proba.initial <- initial.params$startProbs
            size.initial <- initial.params$distributions[,'size']
            prob.initial <- initial.params$distributions[,'prob']
            size.initial[is.na(size.initial)] <- 0
            prob.initial[is.na(prob.initial)] <- 0

        } else if(init == 'random') {
            A.initial <- matrix(stats::runif(numstates^2),ncol=numstates)
            A.initial <- sweep(A.initial,1,rowSums(A.initial),"/")
            proba.initial <- stats::runif(numstates)
            # --> Distributions for dependent states
            size.initial <- stats::runif(1, min=0, max=100) * cumsum(dependent.states.mask)
            # --> Do we assume a mean monosomy between 0 and 100? Isn't this a bit low????
            prob.initial <- stats::runif(1) * dependent.states.mask
            # --> Assign initials for the 0-somy distribution
            index <- which(state.labels == '0-somy')
            size.initial[index] <- 1
            prob.initial[index] <- 0.5

        } else if(init == 'standard') {
            A.initial <- matrix((0.1/(numstates-1)), ncol=numstates, nrow=numstates)
            for(i in 1:numstates)
                A.initial[i,i]  <- 0.9
            proba.initial <- rep(1/numstates, numstates)
            # --> Set initial mean of most.frequent.state distribution to max of count histogram
            max.counts <- as.integer(names(which.max(table(counts[counts>0]))))
            divf <- max(multiplicity[most.frequent.state], 1)
            mean.initial.monosomy <- max.counts/divf
            if(is.na(mean.initial.monosomy))
                mean.initial.monosomy <- 1
            # --> Do you not assume here that the states are consecutive (e.g. 1,2,3 not 1,3,4,...)?
            mean.initial <- mean.initial.monosomy * cumsum(dependent.states.mask)
            var.initial <- mean.initial*2
            size.initial <- rep(0,numstates)
            prob.initial <- rep(0,numstates)
            mask <- dependent.states.mask
            size.initial[mask] <- dnbinom.size(mean.initial[mask], var.initial[mask])
            prob.initial[mask] <- dnbinom.prob(mean.initial[mask], var.initial[mask])
            # --> Assign initials for the 0-somy distribution
            index <- which(state.labels == '0-somy')
            size.initial[index] <- 1
            prob.initial[index] <- 0.5
        }

    hmm <- .C("C_univariate_hmm",
              counts = as.integer(counts),              # --> int* O
              num.bins = as.integer(numbins),           # --> int* T
              num.states = as.integer(numstates),       # --> int* N
              state.labels = as.integer(state.labels),  # --> int* state_labels
              size = double(length=numstates),          # --> double* size
              prob = double(length=numstates),          # --> double* prob
              num.iterations = as.integer(max.iter),    # --> int* maxiter
              time.sec = as.integer(max.time),          # --> double* maxtime
              loglik.delta = as.double(eps.try),        # --> double* eps
              maxPosterior = double(length=numbins),    # --> double* maxPosterior
              states = integer(length=numbins),         # --> int* states
              A = double(length=numstates*numstates),   # --> double* A
              proba = double(length=numstates),         # --> double* proba
              loglik = double(length=1),                # --> double* loglik
              weights = double(length=numstates),       # --> double* weights
              distr.type = as.integer(state.distributions), # --> int* distr_type
              size.initial = as.vector(size.initial),   # --> double* initial_size
              prob.initial = as.vector(prob.initial),   # --> double* initial_prob
              A.initial = as.vector(A.initial),         # --> double* initial_A
              proba.initial = as.vector(proba.initial), # --> double* initial_proba
              use.initial.params = as.logical(1),       # --> bool* use_initial_params
              num.threads = as.integer(num.threads),    # --> int* num_threads
              error = as.integer(0),                    # --> int* error (error handling)
              count.cutoff = as.integer(count.cutoff),  # --> int* count.cutoff
              algorithm = as.integer(algorithm),        # --> int* algorithm
              verbosity = as.integer(verbosity),        # --> int* verbosity
              PACKAGE = 'AneuFinder')

        hmm$eps <- eps.try
        if (num.trials > 1) {
            if(hmm$loglik.delta > hmm$eps)
                warning(paste0("ID = ",ID,": HMM did not converge in trial run ",i_try,"!\n"))
            modellist[[as.character(i_try)]] <- hmm
            init <- 'random'
        } else if(num.trials == 1) {
            if(hmm$loglik.delta > hmm$eps & istep == 1)
                warning(paste0("ID = ",ID,": HMM did not converge!\n"))
        }
    }

    # RERUN MODEL WITH DIFFERENT INITIAL PARAMETERS
    if (num.trials > 1) {
        # --> Mathematically we should select the fit with highest
        # loglikelihood. If we think the fit with the highest loglikelihood is
        # incorrect, we should change the underlying model.
        logliks <- sapply(modellist,'[[','loglik')
        # However, this is very complex and we choose to select a fit that we
        # think is (more) correct, although it has not the highest support
        # given our (imperfect) model.
        df.weight <- as.data.frame(lapply(modellist,'[[','weights'))
        rownames(df.weight) <- state.labels
        # --> Select models where weight of most.frequent.state is at least
        # half of that of actual most frequent state, then select model with
        # highest loglik
        models2use <- (df.weight[most.frequent.state,] / apply(df.weight,2,max)) > 0.5
        models2use[is.na(models2use)] <- FALSE # Needed?

        if (any(models2use))
            index2use <- names(which.max(logliks[models2use]))
        else
            index2use <- names(which.max(logliks))

        hmm <- modellist[[index2use]]
        # --> Check if size and prob parameter are correct
        if (any(is.na(hmm$size) || is.nan(hmm$size) || is.infinite(hmm$size) |
                is.na(hmm$prob) || is.nan(hmm$prob) || is.infinite(hmm$prob)))
            stop("...")

        # --> Rerun the HMM with different epsilon and initial parameters from trial run
        message(paste0("Rerunning trial ",index2use," with eps = ",eps))
        hmm <- .C("C_univariate_hmm",
                  counts = as.integer(counts),             # --> int* O
                  num.bins = as.integer(numbins),          # --> int* T
                  num.states = as.integer(numstates),      # --> int* N
                  state.labels = as.integer(state.labels), # --> int* state_labels
                  size = double(length=numstates),         # --> double* size
                  prob = double(length=numstates),         # --> double* prob
                  num.iterations = as.integer(max.iter),   # --> int* maxiter
                  time.sec = as.integer(max.time),         # --> double* maxtime
                  loglik.delta = as.double(eps),           # --> double* eps
                  maxPosterior = double(length=numbins),   # --> double* maxPosterior
                  states = integer(length=numbins),        # --> int* states
                  A = double(length=numstates*numstates),  # --> double* A
                  proba = double(length=numstates),        # --> double* proba
                  loglik = double(length=1),               # --> double* loglik
                  weights = double(length=numstates),      # --> double* weights
                  distr.type = as.integer(state.distributions), # --> int* distr_type
                  size.initial = as.vector(hmm$size),      # --> double* initial_size
                  prob.initial = as.vector(hmm$prob),      # --> double* initial_prob
                  A.initial = as.vector(hmm$A),            # --> double* initial_A
                  proba.initial = as.vector(hmm$proba),    # --> double* initial_proba
                  use.initial.params = as.logical(1),      # --> bool* use_initial_params
                  num.threads = as.integer(num.threads),   # --> int* num_threads
                  error = as.integer(0),                   # --> int* error (error handling)
                  count.cutoff = as.integer(count.cutoff), # --> int* count.cutoff
                  algorithm = as.integer(algorithm),       # --> int* algorithm
                  verbosity = as.integer(verbosity),       # --> int* verbosity
                  PACKAGE = 'AneuFinder')
    }

    result <- list(ID=ID, bins=binned.data)
    result$bins$state <- state.labels[hmm$states]
    result$bins$copy.number <- multiplicity[as.character(result$bins$state)]
    suppressMessages(result$segments <-
        as(collapseBins(as.data.frame(result$bins),
                        column2collapseBy = 'copy.number',
                        columns2drop = 'width',
                        columns2average = c('counts','mcounts','pcounts')),
           'GRanges'))
    seqlevels(result$segments) <- seqlevels(result$bins) # --> Correct order from as()
    seqlengths(result$segments) <- seqlengths(binned.data)[names(seqlengths(result$segments))]
    result$convergenceInfo <- list(
        eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta,
        num.iterations=hmm$num.iterations ,time.sec=hmm$time.sec, error=hmm$error)
    result$weights <- hmm$weights
    names(result$weights) <- state.labels
    result$startProbs <- hmm$proba
    names(result$startProbs) <- state.labels
    result$startProbs.initial <- hmm$proba.initial
    names(result$startProbs.initial) <- state.labels
    result$transitionProbs <- matrix(hmm$A,ncol=hmm$num.states)
    rownames(result$transitionProbs) <- state.labels
    colnames(result$transitionProbs) <- state.labels
    result$transitionProbs.initial <- matrix(hmm$A.initial,ncol=hmm$num.states)
    rownames(result$transitionProbs.initial) <- state.labels
    colnames(result$transitionProbs.initial) <- state.labels

    distributions <- data.frame()
    distributions.ini <- data.frame()
    for (id in 1:length(hmm$distr.type)) {
        distr <- levels(state.distributions)[hmm$distr.type[id]]

        if (distr == 'dnbinom') {
            distributions <- rbind(distributions, data.frame(
                type = distr,
                size = hmm$size[id],
                prob = hmm$prob[id],
                mu = dnbinom.mean(hmm$size[id], hmm$prob[id]),
                variance = dnbinom.variance(hmm$size[id], hmm$prob[id])))
            distributions.ini <- rbind(distributions.ini, data.frame(
                type = distr,
                size = hmm$size.initial[id],
                prob = hmm$prob.initial[id],
                mu = dnbinom.mean(hmm$size.initial[id],hmm$prob.initial[id]),
                variance = dnbinom.variance(hmm$size.initial[id],hmm$prob.initial[id])))

        } else if(distr == 'dgeom') {
            distributions <- rbind(distributions, data.frame(
                type = distr,
                size = NA,
                prob = hmm$prob[id],
                mu = dgeom.mean(hmm$prob[id]),
                variance = dgeom.variance(hmm$prob[id])))
            distributions.ini <- rbind(distributions.ini, data.frame(
                type=distr,
                size=NA,
                prob=hmm$prob.initial[id],
                mu=dgeom.mean(hmm$prob.initial[id]),
                variance=dgeom.variance(hmm$prob.initial[id])))

        } else if (distr == 'delta') {
            distributions <- rbind(distributions, data.frame(
                type = distr,
                size = NA,
                prob = NA,
                mu = 0,
                variance = 0))
            distributions.ini <- rbind(distributions.ini, data.frame(
                type = distr,
                size = NA,
                prob = NA,
                mu = 0,
                variance = 0))

        } else if (distr == 'dbinom') {
            distributions <- rbind(distributions, data.frame(
                type = distr,
                size = hmm$size[id],
                prob = hmm$prob[id],
                mu = dbinom.mean(hmm$size[id],hmm$prob[id]),
                variance = dbinom.variance(hmm$size[id],hmm$prob[id])))
            distributions.ini <- rbind(distributions.ini, data.frame(
                type = distr,
                size = hmm$size.initial[id],
                prob = hmm$prob.initial[id],
                mu = dbinom.mean(hmm$size.initial[id],hmm$prob.initial[id]),
                variance = dbinom.variance(hmm$size.initial[id], hmm$prob.initial[id])))
        }
    }

    rownames(distributions) <- state.labels
    rownames(distributions.ini) <- state.labels
    result$distributions  <- distributions
    result$distributions.initial <- distributions.ini
    result$warnings <- warlist

    class(result) <- "aneuHMM"
    result
}

initializeStates <- function(states) {
    somy.states <- grep('somy', states, value=TRUE)
    somy.numbers <- as.integer(sapply(strsplit(somy.states, '-somy'), '[[', 1))
    names(somy.numbers) <- somy.states
  
    if ("zero-inflation" %in% states) {
        multiplicity <- c("zero-inflation"=0, somy.numbers)
    } else {
        multiplicity <- somy.numbers
    }
  
    levels.distributions <- c('delta','dgeom','dnbinom','dbinom')
    distributions <- rep(NA, length(states))
    names(distributions) <- states
    distributions[states=='zero-inflation'] <- 'delta'
    distributions[states=='0-somy'] <- 'dgeom'
    distributions[(states != 'zero-inflation') & (states != '0-somy')] <- 'dnbinom'

    states <- factor(states, levels=states)
    distributions <- factor(distributions, levels=levels.distributions)

    list(states=states, distributions=distributions, multiplicity=multiplicity)
}
