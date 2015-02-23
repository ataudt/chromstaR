#' Obtain combinatorial states from specification
#'
#' This function returns all combinatorial (decimal) states that are consistent with a given abstract specification.
#'
#' The binary modification state (unmodified=0 or modified=1) of multiple ChIP-seq samples defines a (decimal) combinatorial state such as:
#' \tabular{ccccccc}{
#'  \tab sample1 \tab sample2 \tab sample3 \tab sample4 \tab sample5 \tab combinatorial state \cr
#' bin1 \tab 0 \tab 0 \tab 1 \tab 0 \tab 0 \tab 4 \cr
#' bin2 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
#' bin3 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \tab 10 \cr
#' bin4 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 15 \cr
#' bin5 \tab 0 \tab 0 \tab 1 \tab 0 \tab 1 \tab 5 \cr
#' }
#'
#' @author Aaron Taudt
#' @param statespec A vector composed of any combination of the following entries: \code{0, 1, 'x', 'r[...]'}, where [\dots] can be any string.
#'   \itemize{
#'     \item \code{0}: sample is 'unmodified'
#'     \item \code{1}: sample is 'modified'
#'     \item \code{'x'}: sample can be both 'unmodified' or 'modified'
#'     \item \code{'rA'}: all samples in group A have to be in the same state
#'     \item \code{'rB'}: all samples in group B have to be in the same state
#'     \item \code{'r[...]'}: all samples in group [\dots] have to be in the same state
#'   }
#' @param diffstatespec A vector composed of any combination of the following entries: \code{'x', 'd[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'x'}: sample can be both 'unmodified' or 'modified'
#'     \item \code{'dA'}: at least one sample in group A has to be different from the other samples in group A 
#'     \item \code{'dB'}: at least one sample in group B has to be different from the other samples in group B 
#'     \item \code{'d[...]'}: at least one sample in group [\dots] has to be different from the other samples in group [\dots] 
#'   }
#' @return A integer vector with (decimal) combinatorial states following the given specification.
#' @examples
#'# Get all combinatorial states where sample1=0, sample2=1, sample3=(0 or 1),
#'#  sample4=sample5
#'stateBrewer(statespec=c(0,1,'x','rA','rA'))
#'
#'# Get all combinatorial states where sample1=sample2=sample3, sample4=sample5
#'stateBrewer(statespec=c('rA','rA','rA','rB','rB'))
#'
#'# Get all combinatorial states where sample1=sample5, sample2=sample3=1,
#'#  sample4=(0 or 1)
#'stateBrewer(statespec=c('rA',1,1,'x','rA'))
#'
#'# Get all combinatorial states where sample1 != sample2, sample3 != sample4,
#'#  sample5=(0 or 1)
#'stateBrewer(statespec=c('x','x','x','x','x'),
#'            diffstatespec=c('dA','dA','dB','dB','x'))
#'# To check whether you chose the correct specification you can use
#'dec2bin(stateBrewer(statespec=c('x','x','x','x','x'),
#'                    diffstatespec=c('dA','dA','dB','dB','x')), ndigits=5)
#' @export
stateBrewer <- function(statespec, diffstatespec=NULL) {

	## Check user input
	for (spec in statespec) {
		if (spec!=1 & spec!=0 & spec!='x' & !grepl('^r', spec)) {
			stop("argument 'statespec' expects a vector composed of any combination of the following entries: 1,0,'x','r[]', where [] can be any string.")
		}
	}
	if (!is.null(diffstatespec)) {
		for (spec in diffstatespec) {
			if (spec!=1 & spec!=0 & spec!='x' & !grepl('^d', spec)) {
				stop("argument 'diffstatespec' expects a vector composed of any combination of the following entries: 'x','d[]', where [] can be any string.")
			}
		}
		if (length(statespec)!=length(diffstatespec)) {
			stop("argument 'diffstatespec' must have the same number of elements as 'statespec'")
		}
	}

	## Variables
	numtracks <- length(statespec)
	groups <- levels(factor(statespec))
	diffgroups <- levels(factor(diffstatespec))

	## Get all possible binary states
	binstates <- dec2bin(0:(2^numtracks-1))

	## Select specified binary states
	for (group in groups) {
		track.index <- which(statespec==group)
		if (group == 0) {
			mask <- !apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('|', x) })
		} else if (group == 1) {
			mask <- apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('&', x) })
		} else if (group == 'x') {
			mask <- rep(T, nrow(binstates))
		} else if (grepl('^r', group)) {
			mask0 <- !apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('|', x) })
			mask1 <- apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('&', x) })
			mask <- mask0 | mask1
		}
		binstates <- binstates[mask,]
	}

	for (diffgroup in diffgroups) {
		track.index <- which(diffstatespec==diffgroup)
		if (grepl('^d', diffgroup)) {
			mask0 <- !apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('|', x) })
			mask1 <- apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('&', x) })
			mask <- !(mask0 | mask1)
		} else if (diffgroup == 'x') {
			mask <- rep(T, nrow(binstates))
		}
		binstates <- binstates[mask,]
	}

	## Convert to decimal
	decstates <- bin2dec(binstates)
	return(decstates)

}
