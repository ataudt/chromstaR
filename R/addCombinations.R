#' Add human interpretable names to combinatorial states
#'
#' This function takes a \code{\link{chromstaR_multivariateHMM}} object and a vector of track names as input and translates the decimal combinatorial states into a human interpretable combination of track names.
#'
#' @author Aaron Taudt
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} object.
#' @param tracknames A vector of track names from which the combinations are generated. This vector must have the same length as the number of tracks in the \code{multi.hmm}. Identical entries will be treated as replicates.
#' @return A \code{\link{chromstaR_multivariateHMM}} object with column 'combinations' added to \code{$bins} and \code{$segments}.
#' @export
addCombinations <- function(multi.hmm, tracknames) {

	if (class(multi.hmm)!=class.multivariate.hmm) {
		stop("argument 'multi.hmm' requires a ", class.multivariate.hmm, " object")
	}
	if (length(tracknames) != length(multi.hmm$IDs)) {
		stop("argument 'tracknames' needs to have length(multi.hmm$IDs) = ", length(multi.hmm$IDs), " elements")
	}
	statespec <- paste0('r.', tracknames)
	mapping <- stateBrewer(statespec, inverse=TRUE)
	multi.hmm$bins$combination <- factor(mapping[as.character(multi.hmm$bins$state)])
	multi.hmm$segments$combination <- factor(mapping[as.character(multi.hmm$segments$state)])

	if (any(is.na(multi.hmm$segments$combination))) {
		stop("NAs produced in combinations. You probably made a mistake with argument 'tracknames'. Keep in mind that identical entries will be treated as replicates with identical peak calls.")
	}

	return(multi.hmm)

}
