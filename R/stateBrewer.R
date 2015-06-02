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
#' @param statespec A vector composed of any combination of the following entries: \code{'0.[]', '1.[]', 'x.[]', 'r.[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'0.A'}: sample A is 'unmodified'
#'     \item \code{'1.B'}: sample B is 'modified'
#'     \item \code{'x.C'}: sample C can be both 'unmodified' or 'modified'
#'     \item \code{'r.D'}: all samples in group D have to be in the same state
#'     \item \code{'r.[]'}: all samples in group [] have to be in the same state
#'   }
#' @param diffstatespec A vector composed of any combination of the following entries: \code{'x.[]', 'd.[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'x.A'}: sample A can be both 'unmodified' or 'modified'
#'     \item \code{'d.B'}: at least one sample in group B has to be different from the other samples in group A 
#'     \item \code{'d[]'}: at least one sample in group [] has to be different from the other samples in group [] 
#'   }
#' @param inverse If \code{TRUE}, names and entries of the output are swapped.
#' @param sep Separator used to separate the tracknames in the combinations.
#' @return A named integer vector with (decimal) combinatorial states following the given specification. The names of the vector are the combinatorial states composed of the different groups. If \code{inverse=TRUE}, names and numbers are swapped.
#' @examples
#'# Get all combinatorial states where sample1=0, sample2=1, sample3=(0 or 1),
#'#  sample4=sample5
#'stateBrewer(statespec=c('0.A','1.B','x.C','r.D','r.D'))
#'
#'# Get all combinatorial states where sample1=sample2=sample3, sample4=sample5
#'stateBrewer(statespec=c('r.A','r.A','r.A','r.B','r.B'))
#'
#'# Get all combinatorial states where sample1=sample5, sample2=sample3=1,
#'#  sample4=(0 or 1)
#'stateBrewer(statespec=c('r.A','1.B','1.C','x.D','r.A'))
#'
#'# Get all combinatorial states where sample1 != sample2, sample3 != sample4,
#'#  sample5=(0 or 1)
#'stateBrewer(statespec=c('x.A','x.B','x.C','x.D','x.E'),
#'            diffstatespec=c('d.a','d.a','d.b','d.b','x.c'))
#'# To check whether you chose the correct specification you can use
#'dec2bin(stateBrewer(statespec=c('x.A','x.B','x.C','x.D','x.E'),
#'            diffstatespec=c('d.a','d.a','d.b','d.b','x.c')), colnames=c('A','B','C','D','E'))
#' @export
stateBrewer <- function(statespec, diffstatespec=NULL, inverse=FALSE, sep='-') {

	## Check user input
	for (spec in statespec) {
		if (!grepl('^1\\.', spec) & !grepl('^0\\.', spec) & !grepl('^x\\.', spec) & !grepl('^r\\.', spec)) {
			stop("argument 'statespec' expects a vector composed of any combination of the following entries: '1.[]','0.[]','x.[]','r.[]', where [] can be any string.")
		}
	}
	if (!is.null(diffstatespec)) {
		for (spec in diffstatespec) {
			if (!grepl('^1\\.', spec) & !grepl('^0\\.', spec) & !grepl('^x\\.', spec) & !grepl('^d\\.', spec)) {
				stop("argument 'diffstatespec' expects a vector composed of any combination of the following entries: '1.[]','0.[]','x.[]','r.[]', where [] can be any string.")
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
	tracknames <- sub('^.\\.', '', statespec)

	## Generate specified binary states
	numstates <- 2^length(which(!grepl('^0\\.|^1\\.', groups)))
	binstates <- matrix(FALSE, ncol=numtracks, nrow=numstates)
	i1 <- 1
	for (group in groups) {
		track.index <- which(statespec==group)
		for (itrack in track.index) {
			if (grepl('^1\\.', group)) {
				binstates[,itrack] <- TRUE
			} else if (grepl('^x\\.', group)) {
				numeach <- numstates/2 / 2^(i1-1)
				binstates[,itrack] <- rep(c(rep(FALSE, numeach), rep(TRUE, numeach)), 2^(i1-1))
				i1 <- i1 + 1
			} else if (grepl('^r\\.', group)) {
				numeach <- numstates/2 / 2^(i1-1)
				binstates[,itrack] <- rep(c(rep(FALSE, numeach), rep(TRUE, numeach)), 2^(i1-1))
			}
		}
		if (grepl('^r\\.', group)) {
			i1 <- i1 + 1
		}
	}
	colnames(binstates) <- tracknames

	for (diffgroup in diffgroups) {
		track.index <- which(diffstatespec==diffgroup)
		if (grepl('^d\\.', diffgroup)) {
			mask0 <- !apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('|', x) })
			mask1 <- apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('&', x) })
			mask <- !(mask0 | mask1)
		} else if (grepl('^x\\.', diffgroup)) {
			mask <- rep(T, nrow(binstates))
		}
		binstates <- binstates[mask,]
		if (class(binstates)!='matrix') {
			binstates <- matrix(binstates, ncol=length(binstates))
			colnames(binstates) <- tracknames
		}
	}

	## Construct state names
	mask <- !grepl('^r\\.', statespec) | !duplicated(statespec) # only first replicate
	tracknames.mask <- tracknames[mask]
	if (nrow(binstates) > 1) {
		statenames <- unlist(lapply(apply(as.matrix(binstates[,mask]), 1, function(x) { tracknames.mask[x] }), paste, collapse=sep))
	} else {
		statenames <- paste(tracknames.mask[binstates[,mask]], collapse=sep)
	}
	## Convert to decimal
	decstates <- bin2dec(binstates)
	names(decstates) <- statenames

	## Change name and entries
	if (inverse) {
		decstates.inverse <- names(decstates)
		names(decstates.inverse) <- decstates
		decstates <- decstates.inverse
	}

	return(decstates)

}
