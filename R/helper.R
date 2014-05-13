# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- c("zero-inflation","unmodified","modified")
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.chromstar.univariate <- "chromstar.univariate.hmm"
class.chromstar.multivariate <- "chromstar.multivariate.hmm"
 
# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
fsize <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

fprob <- function(mean, variance) {
	return(mean/variance)
}

fmean <- function(size, prob) {
	return(size/prob - size)
}

fvariance <- function(size, prob) {
	return( (size - prob*size) / prob^2 )
}
