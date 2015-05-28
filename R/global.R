#' @useDynLib chromstaR
#' @import GenomicRanges
#' @import IRanges
NULL

# =======================================================
# Some global variables that can be used in all functions
# =======================================================
state.labels <- factor(c("zero-inflation","unmodified","modified"), levels=c("zero-inflation","unmodified","modified"))
state.distributions <- factor(c('delta','dnbinom','dnbinom'), levels=c('delta','dnbinom'))
class.univariate.hmm <- "chromstaR_univariateHMM"
class.multivariate.hmm <- "chromstaR_multivariateHMM"
state.colors <- c("zero-inflation"="gray30","unmodified"="gray48","modified"="orangered3", "total"="black", "reads"="grey35")

get.state.labels <- function() { return(state.labels) }

#' Get state colors
#'
#' Get the colors that are used for plotting
#' @export
get.state.colors <- function() { return(state.colors[state.labels]) }
 
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

