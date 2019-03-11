#' Find copy number variations (bivariate)
#'
#' \code{biHMM.findCNVs} finds CNVs using read count information from both
#' strands.
#'
#' @inheritParams HMM.findCNVs
#' @inheritParams findCNVs
#' @return An \code{\link{aneuBiHMM}} object.
#' @importFrom stats pgeom pnbinom qnorm
biHMM.findCNVs <- function(binned.data, ID=NULL, eps=0.01, init="standard",
                           max.time=-1, max.iter=-1, num.trials=1,
                           eps.try=NULL, num.threads=1,
                           count.cutoff.quantile=0.999,
                           states=c("zero-inflation",paste0(0:10,"-somy")),
                           most.frequent.state="1-somy", algorithm='EM',
                           initial.params=NULL, verbosity=1) {

    on.exit(.C("C_multivariate_cleanup", as.integer(num.comb.states),
               PACKAGE='AneuFinder'))

    warlist <- list()
    if (!is.null(initial.params))
        init <- 'initial.params'
    if (is.null(eps.try))
        eps.try <- eps

    counts <- matrix(
        c(mcols(binned.data)[,'mcounts'], mcols(binned.data)[,'pcounts']), ncol=2,
        dimnames=list(bin=1:length(binned.data), strand=c('minus','plus')))
    count.cutoff <- ceiling(quantile(counts,count.cutoff.quantile))
    rows.higher <- which(counts > count.cutoff)
    counts[rows.higher] <- count.cutoff

    if (length(rows.higher) > 0)
        message("Replaced read counts > ",count.cutoff," (quantile",
                count.cutoff.quantile, ") by ", count.cutoff, " in ",
                length(rows.higher), " bins. Set option 'count.cutoff.quantile=1' to",
                "disable this filtering.\n", "This filtering was done to enhance
                performance.")

    if (any(is.na(counts))) {
        stop(paste0("ID = ",ID,": NAs found in reads."))
    } else if(all(counts == 0)) {
        wstr = paste0("ID = ",ID,": All counts in data are zero. No HMM done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings <- warlist
        return(result)
    } else if (any(counts < 0)) {
        wstr = paste0("ID = ",ID,": Some counts in data are negative. No HMM done.")
        warlist[[length(warlist)+1]] <- warning(wstr)
        result$warnings <- warlist
        return(result)
    }

    if (init == 'initial.params') {
        uni.transitionProbs <- initial.params$univariateParams$transitionProbs
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

    } else if(init == 'standard') {
        binned.data.minus   <- binned.data # --> Stack the strands and run HMM.findCNVs
        strand(binned.data.minus) <- '-'
        binned.data.minus$counts <- binned.data.minus$mcounts
        binned.data.plus    <- binned.data
        strand(binned.data.plus) <- '+'
        binned.data.plus$counts <- binned.data.plus$pcounts
        binned.data.stacked <- c(binned.data.minus,binned.data.plus)

        message("Running univariate HMM")
        # --> We here run an HMM on two sets of bins after eachother. A bit weird?????????
        model.stacked <- HMM.findCNVs(binned.data.stacked, ID, eps=eps,
                                      init=init, max.time=max.time,
                                      max.iter=max.iter, num.trials=num.trials,
                                      eps.try=eps.try, num.threads=num.threads,
                                      count.cutoff.quantile=1, states=states,
                                      most.frequent.state=most.frequent.state)
        if (is.na(model.stacked$convergenceInfo$error)) {
            result$warnings <- model.stacked$warnings # --> What if not run? still an NA?? I think NULL
            return(result)
        }
        uni.transitionProbs <- model.stacked$transitionProbs # --> Extract probs, distributions
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
        comb.states.per.bin <-
            paste(model.stacked$bins$state[as.character(strand(model.stacked$bins))=='-'],
                  model.stacked$bins$state[as.character(strand(model.stacked$bins))=='+'])
        comb.states.per.bin <- factor(comb.states.per.bin, levels=comb.states)
        A.initial           <- double(length=num.comb.states^2)
        proba.initial       <- double(length=num.comb.states)
        use.initial         <- FALSE
    }

    # CALCULATE DENSITIES THAT ARE NEEDED AS INPUT FOR MULTIVARIATE HMM
    # --> Pre-compute z-values for each number of counts
    z.per.count <- array(NA, dim=c(count.cutoff+1,num.uni.states),
                         dimnames=list(counts=0:count.cutoff,state=uni.states))
    xcounts <- 0:count.cutoff
    for (istate in 1:num.uni.states) {
        if (distributions[istate,'type'] == 'dnbinom') {
            size  <- distributions[istate,'size']
            prob  <- distributions[istate,'prob']
            u     <- stats::pnbinom(xcounts,size,prob)
        } else if (distributions[istate,'type'] == 'dbinom') {
            size  <- distributions[istate,'size']
            prob  <- distributions[istate,'prob']
            u     <- stats::pbinom(xcounts,size,prob)
        } else if (distributions[istate,'type'] == 'delta') {
            u     <- rep(1, length(xcounts))
        } else if (distributions[istate,'type'] == 'dgeom') {
            prob  <- distributions[istate,'prob']
            u     <- stats::pgeom(xcounts, prob)
        }
        qnorm_u <- stats::qnorm(u)
        mask <- qnorm_u==Inf # --> What about -Inf??????
        qnorm_u[mask] <- stats::qnorm(1-1e-16) # --> Why???
        z.per.count[,istate] <- qnorm_u
    }

    num.bins <- length(binned.data)
    z.per.bin <- array(NA, dim=c(num.bins,num.models,num.uni.states),
                       dimnames=list(bin=1:num.bins,strand=c("minus","plus"),
                                     state=uni.states))
    for (istrand in 1:num.models)
        for (istate in 1:num.uni.states)
            z.per.bin[,istrand,istate] <- z.per.count[counts[,istrand]+1,istate]

    correlationMatrix <- array(0, dim=c(num.models,num.models,num.comb.states),
                               dimnames=list(strand=c("minus","plus"),
                                             strand=c("minus","plus"),
                                             comb.state=comb.states))
    correlationMatrixInverse <- correlationMatrix
    det.val <- rep(0,num.comb.states)
    names(det.val) <- comb.states

    for (comb.state in comb.states) {
        state <- strsplit(comb.state,' ')[[1]]
        mask <- which(comb.states.per.bin == comb.state)
        z.temp <- matrix(NA,ncol=num.models,nrow=length(mask))

        for(istrand in 1:num.models)
            z.temp[,istrand]  <- z.per.bin[mask,istrand,state[istrand]]

        temp <- tryCatch({
            if (nrow(z.temp) > 1) {
                correlationMatrix[,,comb.state] <- cor(z.temp)
                det.val[comb.state] <- det(correlationMatrix[,,comb.state])
                correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
            } else {
                correlationMatrix[,,comb.state] <- diag(num.models)
                det.val[comb.state] <- det(correlationMatrix[,,comb.state])
                correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
            }
            0
        }, warning = function(war) 1, error = function(err) 0)

        if (temp != 0) { # error or warning
            correlationMatrix[,,comb.state] <- diag(num.models)
            det.val[comb.state] <- det(correlationMatrix[,,comb.state])
            correlationMatrixInverse[,,comb.state] <- solve(correlationMatrix[,,comb.state])
        }
    }

    # --> Calculate multivariate densities for each state
    # --> At the end density has very small values!!!!! even e-300!!!!! Does this ever go wrong???
    densities <- matrix(1, ncol=num.comb.states, nrow=num.bins,
                        dimnames=list(bin=1:num.bins,comb.state=comb.states))

    for (comb.state in comb.states) {
        istate              <- which(comb.state==comb.states)
        state               <- strsplit(comb.state,' ')[[1]]
        z.temp              <- matrix(NA,ncol=num.models,nrow=num.bins)
        product             <- 1
        for (istrand in 1:num.models){
            z.temp[,istrand] <- z.per.bin[,istrand,state[istrand]]
            if (distributions[state[istrand],'type'] == 'dnbinom') {
                size         <- distributions[state[istrand],'size']
                prob         <- distributions[state[istrand],'prob']
                product      <- product * stats::dnbinom(counts[,istrand],size,prob)
            } else if(distributions[state[istrand],'type'] == 'dgeom') {
                prob         <- distributions[state[istrand],'prob']
                product      <- product * stats::dgeom(counts[,istrand],prob)
            } else if(distributions[state[istrand],'type'] == 'delta') {
                product      <- product * ifelse(counts[,istrand]==0, 1, 0)
            }
        }
        exponent <- -0.5 * apply((z.temp %*% (correlationMatrixInverse[,,istate] - diag(num.models)))*z.temp,
                             1, sum)
        exponent[is.nan(exponent)] <- 0
        densities[,istate]  <- product * det.val[istate]^(-0.5) * exp(exponent)
    }

    densities[densities > 1] <- 1 # --> Check if densities are > 1 or < 0
    densities[densities < 0] <- 0 # --> Is it ok if this happens???????
    check <- which(apply(densities, 1, sum) == 0) # --> Check if densities are 0 everywhere in some bins

    if (length(check) > 0) {
        if(check[1] == 1){
            densities[1,] <- rep(1e-10, ncol(densities)) # --> Why 10**-10??????
            check <- check[-1]
        }
        # --> If row sum is zero then take values of previous row. Why?????????
        for (icheck in check)
            densities[icheck,] <- densities[icheck-1,]
    }

    # RUN THE MULTIVARIATE HMM
    algorithm <- factor(algorithm, levels=c('baumWelch','viterbi','EM'))
    hmm <- .C("C_multivariate_hmm",
              densities = as.double(densities),         # --> double* D
              num.bins = as.integer(num.bins),          # --> int* T
              num.comb.states = as.integer(num.comb.states), # --> int* N
              num.strands = as.integer(num.models),     # --> int* Nmod
              comb.states = as.integer(comb.states),    # --> int* comb_states
              num.iterations = as.integer(max.iter),    # --> int* maxiter
              time.sec = as.integer(max.time),          # --> double* maxtime
              loglik.delta = as.double(eps),            # --> double* eps
              maxPosterior = double(length=num.bins),   # --> double* maxPosterior
              states = integer(length=num.bins),        # --> int* states
              A = double(length=num.comb.states*num.comb.states), # --> double* A
              proba = double(length=num.comb.states),   # --> double* proba
              loglik = double(length=1),                # --> double* loglik
              A.initial = as.vector(A.initial),         # --> double* initial_A
              proba.initial = as.vector(proba.initial), # --> double* initial_proba
              use.initial.params = as.logical(use.initial), # --> bool* use_initial_params
              num.threads = as.integer(num.threads),    # --> int* num_threads
              error = as.integer(0),                    # --> error handling
              algorithm = as.integer(algorithm),        # --> int* algorithm
              verbosity = as.integer(verbosity),        # --> int* verbosity
              PACKAGE = 'AneuFinder')

    if(hmm$loglik.delta > eps)
        warning(paste0("ID = ",ID,": HMM did not converge!\n"))

    if(hmm$error == 1) {
        warlist[[length(warlist)+1]] <-
            warning("ID = ", ID, ": A NaN occurred during the Baum-Welch! ",
                    "Parameter estimation terminated prematurely.\n",
                    "Check your library! The following factors are known to cause this error:\n",
                    "1) Your read counts contain very high numbers.",
                    " Try again with a lower value for 'count.cutoff.quantile'.\n",
                    "2) Your library contains too few reads in each bin.\n",
                    "3) Your library contains reads for a different genome than it was aligned to.")
        result$warnings <- warlist
        return(result)
    } else if(hmm$error == 2) {
        warlist[[length(warlist)+1]] <- warning("ID = ", ID,
            ": An error occurred during the Baum-Welch!\n",
            "Parameter estimation terminated prematurely. Check your library.")
        result$warnings     <- warlist
        return(result)
    }

    result                <- list(ID=ID, bins=binned.data)
    matrix.states         <- do.call("rbind",strsplit(as.character(comb.states[hmm$states]),split=' '))
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
    result$bins$combi     <- paste(result$bins$mcopy.number, result$bins$pcopy.number)
    suppressMessages(result$segments <-
        as(collapseBins(as.data.frame(result$bins),
                        column2collapseBy = 'combi',
                        columns2drop = 'width',
                        columns2average = c('counts','mcounts','pcounts')),
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
    result$bins$combi <- NULL
    result$convergenceInfo <- list(eps=eps,loglik=hmm$loglik,loglik.delta=hmm$loglik.delta,
                                 num.iterations=hmm$num.iterations,time.sec=hmm$time.sec)
    result$weights <- table(result$bins$state)/length(result$bins)
    result$startProbs <- hmm$proba
    names(result$startProbs) <- comb.states
    result$startProbs.initial <- hmm$proba.initial
    names(result$startProbs.initial) <- comb.states
    result$transitionProbs <- matrix(hmm$A, ncol=num.comb.states)
    colnames(result$transitionProbs) <- comb.states
    rownames(result$transitionProbs) <- comb.states
    result$transitionProbs.initial <- matrix(hmm$A.initial, ncol=num.comb.states)
    colnames(result$transitionProbs.initial) <- comb.states
    rownames(result$transitionProbs.initial) <- comb.states
    result$distributions  <- distributions
    result$univariateParams <- list(weights=uni.weights,
                                    startProbs=uni.startProbs,
                                    transitionProbs=uni.transitionProbs,
                                    distributions=distributions)
    result$warnings <- warlist

    class(result) <- "aneuBiHMM"
    result
}
