# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)

dnbinom.size <- function(mean, variance) {
	mean^2 / (variance - mean)
}

dnbinom.prob <- function(mean, variance) {
	mean/variance
}

dnbinom.mean <- function(size, prob) {
	size/prob - size
}

dnbinom.variance <- function(size, prob) {
	(size - prob*size) / prob^2
}

dgeom.prob <- function(mean) {
  1/(1+mean)
}

dgeom.mean <- function(prob) {
	(1-prob)/prob
}

dgeom.variance <- function(prob) {
	(1-prob)/prob^2
}

dbinom.size <- function(mean, variance) {
	mean^2/(mean-variance)
}

dbinom.prob <- function(mean, variance) {
	(mean-variance)/mean
}

dbinom.mean <- function(size, prob) {
	size*prob
}

dbinom.variance <- function(size, prob) {
	size*prob * (1-prob)
}

assign.distributions <- function(counts, states) {
    distributions <- list()
    # --> Determine distribution type based on mu and variance.
    bins.splt <- split(counts,states)

    # --> Is this correct?????? Why do we need this? Originally the results
    # were not based on these distributions, right?
    for (i1 in 1:length(bins.splt)) {
        if (length(bins.splt[[i1]]) > 0) {
            qus <- quantile(bins.splt[[i1]], c(0.01, 0.99))
            qcounts <- bins.splt[[i1]]
            qcounts <- qcounts[qcounts >= qus[1] & qcounts <= qus[2]]
            if(sum(qcounts) == 0 | length(qcounts)==0)
                qcounts <- bins.splt[[i1]]

            mu <- mean(qcounts)

            variance <- var(qcounts)
            if (is.na(variance))
                variance <- mu + 1 # --> Somewhat arbitrary. A bit strange?????

            if (names(bins.splt)[i1] == '0-somy') {
                distr <- 'dgeom'
                size <- NA
                prob <- dgeom.prob(mu)
            } else {
                if (is.na(variance) | is.na(mu)) {
                    distr <- 'dnbinom'
                    size <- NA
                    prob <- NA
                } else {
                    if (variance < mu) {
                        distr <- 'dbinom'
                        size <- dbinom.size(mu,variance)
                        prob <- dbinom.prob(mu,variance)
                    } else if(variance > mu) {
                        distr <- 'dnbinom'
                        size <- dnbinom.size(mu,variance)
                        prob <- dnbinom.prob(mu,variance)
                    } else {
                        distr <- 'dpois'
                        size <- NA
                        prob <- mu
                    }
                }
            }
            distributions[[i1]] <- data.frame(type=distr, size=size, prob=prob,
                                              mu=mu, variance=variance)
        }
    }

    distributions <- do.call(rbind, distributions)
    rownames(distributions) <- names(which(lapply(bins.splt,length) > 0))
    distributions
}
