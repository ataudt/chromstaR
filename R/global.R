#' @useDynLib chromstaR, .registration = TRUE, .fixes = ""
#' @import GenomeInfoDb
#' @import IRanges
#' @import GenomicRanges
NULL

# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- factor(c("zero-inflation","unmodified","modified"), levels=c("zero-inflation","unmodified","modified"))
state.distributions <- factor(c('delta','dnbinom','dnbinom'), levels=c('delta','dnbinom'))
class.univariate.hmm <- "uniHMM"
class.multivariate.hmm <- "multiHMM"
class.combined.multivariate.hmm <- "combinedMultiHMM"

# ============================================================================
# Functions for a Negative Binomial to transform (mean,variance)<->(size,prob)
# ============================================================================
dnbinom.size <- function(mean, variance) {
	return(mean^2 / (variance - mean))
}

dnbinom.prob <- function(mean, variance) {
	return(mean/variance)
}

dnbinom.mean <- function(size, prob) {
	return(size/prob - size)
}

dnbinom.variance <- function(size, prob) {
	return( (size - prob*size) / prob^2 )
}

# ==================
# Printing functions
# ==================
message.underlined <- function(string, line='=') {
	string.length <- nchar(string)
	underline <- paste0(rep(line, string.length), collapse='')
	message(string)
	message(underline)
}

message.overlined <- function(string, line='=') {
	string.length <- nchar(string)
	overline <- paste0(rep(line, string.length), collapse='')
	message(string)
	message(overline)
}

