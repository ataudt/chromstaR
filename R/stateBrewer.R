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
#' @param diff.conditions A vector with the same length as \code{statespec}. Similar entries will be treated as belonging to the same condition. If this parameter is specified, only states that are different between the conditions are returned.
#' @param tracks2compare A vector with the same length as \code{statespec}. This vector defines the tracks between which conditions are compared.
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
#' @export
stateBrewer <- function(statespec, diff.conditions=NULL, tracks2compare=NULL, mindiff=1, inverse=FALSE, sep='-') {

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
	tracknames <- sub('^.\\.', '', statespec)

	### Generate specified binary states ###
	groups <- levels(factor(statespec))
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

	### Select differential states ###
	if (!is.null(diff.conditions)) {
		if (is.null(tracks2compare)) {
			stop("argument 'tracks2compare' must be specified if 'diff.conditions' was specified")
		}
		diffgroups <- unique(diff.conditions)
		tracks2compare.split <- split(tracks2compare, diff.conditions)
		intersect.tracks <- Reduce(intersect, lapply(tracks2compare.split, unique))
		bindiffmatrix <- dec2bin(0:(2^length(intersect.tracks)-1))
		controlsum <- apply(bindiffmatrix, 1, sum)
		bindiffmatrix <- bindiffmatrix[controlsum >= mindiff,]
		if (class(bindiffmatrix)!='matrix') {
			bindiffmatrix <- matrix(bindiffmatrix, nrow=1)
		}
		diffstatespec.list <- list()
		for (tracks in tracks2compare.split) {
			#TODO: tracksNOT2use
			tracks2use <- tracks[tracks %in% intersect.tracks]
			num.tracks.split <- length(tracks2use)
			num.reps <- rle(as.integer(factor(tracks, levels=unique(tracks2use))))$lengths
			cum.num.reps <- cumsum(num.reps)
			diffstatespec.part <- apply(bindiffmatrix, 1, function(x) { c('x.','d.')[x+1] })
			diffstatespec.part.reps <- matrix(NA, ncol=sum(num.reps), nrow=nrow(bindiffmatrix))
			num.rep_prev <- 1
			for (i1 in 1:length(cum.num.reps)) {
				diffstatespec.part.reps[,num.rep_prev:cum.num.reps[i1]] <- rep(diffstatespec.part[i1,], num.reps[i1])
				num.rep_prev <- cum.num.reps[i1]+1
			}
			diffstatespec.list[[length(diffstatespec.list)+1]] <- t(apply(diffstatespec.part.reps, 1, function(x) { paste0(x, tracks2use) }))
		}
		diffstatespecs <- do.call(cbind, diffstatespec.list)

		## Go through all diffstate specifications
		binstates.list <- list()
		for (irow in 1:nrow(diffstatespecs)) {
			diffstatespec <- diffstatespecs[irow,]
			diffgroups <- levels(factor(diffstatespec))
			binstates.irow <- binstates
			for (diffgroup in diffgroups) {
				track.index <- which(diffstatespec==diffgroup)
				if (grepl('^d\\.', diffgroup)) {
					mask0 <- !apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('|', x) })
					mask1 <- apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('&', x) })
					mask <- !(mask0 | mask1)
				} else if (grepl('^x\\.', diffgroup)) {
					mask <- rep(T, nrow(binstates.irow))
				}
				binstates.irow <- binstates.irow[mask,]
				if (class(binstates.irow)!='matrix') {
					binstates.irow <- matrix(binstates.irow, ncol=length(binstates.irow))
					colnames(binstates.irow) <- tracknames
				}
			}
			binstates.list[[irow]] <- binstates.irow
		}
		binstates <- do.call(rbind, binstates.list)
	}

	### Construct state names ###
	mask <- !grepl('^r\\.', statespec) | !duplicated(statespec) # only first replicate
	tracknames.mask <- tracknames[mask]
	if (nrow(binstates) > 1) {
		statenames.sep <- apply(as.matrix(binstates[,mask]), 1, function(x) { tracknames.mask[x] })
		if (class(statenames.sep)=='list') {
			statenames <- unlist(lapply(statenames.sep, paste, collapse=sep))
		} else if (class(statenames.sep)=='matrix') {
			statenames <- apply(statenames.sep, 2, paste, collapse=sep)
		}
	} else {
		statenames <- paste(tracknames.mask[binstates[,mask]], collapse=sep)
	}
	## Convert to decimal
	decstates <- bin2dec(binstates)
	duplicate.mask <- !duplicated(decstates)
	decstates <- decstates[duplicate.mask]
	names(decstates) <- statenames[duplicate.mask]

	## Change name and entries
	if (inverse) {
		decstates.inverse <- names(decstates)
		names(decstates.inverse) <- decstates
		decstates <- decstates.inverse
	}

	return(decstates)

}
