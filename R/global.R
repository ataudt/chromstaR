# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- c("zero-inflation","unmodified","modified")
state.distributions <- c('delta','dnbinom','dnbinom')
coordinate.names <- c("chrom","start","end")
binned.data.names <- c(coordinate.names,"reads")
class.chromstar.univariate <- "chromstar.univariate.hmm"
class.chromstar.multivariate <- "chromstar.multivariate.hmm"
gcolors <- c("zero-inflation"="gray30","unmodified"="gray48","modified"="orangered3", "total"="black")
 
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

# ====================================================
# Estimate the RAM consumption of the multivariate HMM
# ====================================================
estimate.RAM.consumption <- function(num.states, num.bins, num.mods) {
# scalefactoralpha T
# scalealpha N,T
# scalebeta N,T
# densities N,T
# gamma N,T
# A N,N
# sumxi N,N
# gammaold N,T
# multiO Nmod,T int
# binary_states N,Nmod bool

# double 8bytes
# int 4bytes
# bool 4bytes
	ram <- 8 * (num.bins + num.bins*num.states *5 + num.states^2 *2) + 4 * (num.mods*num.bins) + 4 * (num.states*num.mods)
	ram.Mb <- as.integer(ram / 2^20)
	return(ram.Mb)
}
