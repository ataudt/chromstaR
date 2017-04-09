#' @useDynLib chromstaR, .registration = TRUE, .fixes = ""
#' @import GenomeInfoDb
#' @import IRanges
#' @import GenomicRanges
#' @import chromstaRData
#' @importFrom methods as is
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

rpkm.vector <- function(counts, binsize) {
    rpkm <- counts / sum(as.numeric(counts)) * 1e6 * 1000 / binsize
    rpkm[is.nan(rpkm)] <- 0
    return(rpkm)
}

rpkm.matrix <- function(counts, binsize) {
    rpkm <- sweep(counts, MARGIN = 2, STATS = colSums(counts), FUN = '/')
    rpkm[is.nan(rpkm)] <- 0
    rpkm <- rpkm * 1e6 * 1000 / binsize
    return(rpkm)
}