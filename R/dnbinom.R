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
