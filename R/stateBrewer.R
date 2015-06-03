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
#' @param replicates A vector specifying the replicate structure. Similar entries will be treated as replicates.
#' @param differential.states A logical specifying whether differential states shall be returned.
#' @param min.diff The minimum number of differences between conditions.
#' @param common.states A logical specifying whether common states shall be returned.
#' @param conditions A vector with the same length as \code{replicates}. Similar entries will be treated as belonging to the same condition. If this parameter is specified, only states that are different between the conditions are returned.
#' @param tracks2compare A vector with the same length as \code{replicates}. This vector defines the tracks between which conditions are compared.
#' @param inverse If \code{TRUE}, names and entries of the output are swapped.
#' @param sep Separator used to separate the tracknames in the combinations.
#' @param statespec If this parameter is specified, \code{replicates} will be ignored. A vector composed of any combination of the following entries: \code{'0.[]', '1.[]', 'x.[]', 'r.[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'0.A'}: sample A is 'unmodified'
#'     \item \code{'1.B'}: sample B is 'modified'
#'     \item \code{'x.C'}: sample C can be both 'unmodified' or 'modified'
#'     \item \code{'r.D'}: all samples in group D have to be in the same state
#'     \item \code{'r.[]'}: all samples in group [] have to be in the same state
#'   }
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
stateBrewer <- function(replicates=NULL, differential.states=FALSE, min.diff=1, common.states=FALSE, conditions=NULL, tracks2compare=NULL, inverse=FALSE, sep='-', statespec=NULL) {

	## Check user input
	if (is.null(statespec)) {
		if (!is.null(replicates)) {
			statespec <- paste0('r.', replicates)
		} else {
			stop("Please specify either 'replicates' or 'statespec'.")
		}
	}
	for (spec in statespec) {
		if (!grepl('^1\\.', spec) & !grepl('^0\\.', spec) & !grepl('^x\\.', spec) & !grepl('^r\\.', spec)) {
			stop("argument 'statespec' expects a vector composed of any combination of the following entries: '1.[]','0.[]','x.[]','r.[]', where [] can be any string.")
		}
	}
	if (differential.states | common.states) {
		if (is.null(conditions) | is.null(tracks2compare)) {
			stop("Please specify 'conditions' and 'tracks2compare' if you want to obtain differential or common states.")
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
	binstates.diff <- NULL
	if (differential.states) {
		if (is.null(tracks2compare)) {
			stop("argument 'tracks2compare' must be specified if 'conditions' was specified")
		}
		tracks2compare.split <- split(tracks2compare, conditions)
		intersect.tracks <- Reduce(intersect, lapply(tracks2compare.split, unique))
		bindiffmatrix <- dec2bin(0:(2^length(intersect.tracks)-1))
		controlsum <- apply(bindiffmatrix, 1, sum)
		bindiffmatrix <- bindiffmatrix[controlsum >= min.diff,]
		if (class(bindiffmatrix)!='matrix') {
			bindiffmatrix <- matrix(bindiffmatrix, nrow=1)
		}
		diffstatespec.list <- list()
		for (tracks in tracks2compare.split) {
			#TODO: tracksNOT2use
			tracks2use <- tracks[tracks %in% intersect.tracks]
			num.tracks.split <- length(tracks2use)
			diffstatespec.part <- t(apply(bindiffmatrix, 1, function(x) { c('x.','d.')[x+1] }))
			colnames(diffstatespec.part) <- intersect.tracks
			diffstatespec.part.reps <- matrix(NA, ncol=length(tracks2use), nrow=nrow(bindiffmatrix))
			for (track in intersect.tracks) {
				index <- which(track==tracks2use)
				diffstatespec.part.reps[,index] <- rep(diffstatespec.part[,as.character(track)], length(index))
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
					mask0 <- !apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('|', x) }) # rows where all group members are 0
					mask1 <- apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('&', x) }) # rows where all group members are 1
					mask <- !(mask0 | mask1) # rows where not all group members are either 0 or 1
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
		binstates.diff <- do.call(rbind, binstates.list)
	}
	# There are still duplicate rows at this point, they are removed at the end

	### Select common states ###
	binstates.common <- NULL
	if (common.states) {
		if (is.null(tracks2compare)) {
			stop("argument 'tracks2compare' must be specified if 'conditions' was specified")
		}
		tracks2compare.split <- split(tracks2compare, conditions)
		intersect.tracks <- Reduce(intersect, lapply(tracks2compare.split, unique))
		bincommonmatrix <- dec2bin(0:(2^length(intersect.tracks)-1))
		if (class(bincommonmatrix)!='matrix') {
			bincommonmatrix <- matrix(bincommonmatrix, nrow=1)
		}
		commonstatespec.list <- list()
		for (tracks in tracks2compare.split) {
			#TODO: tracksNOT2use
			tracks2use <- tracks[tracks %in% intersect.tracks]
			num.tracks.split <- length(tracks2use)
			commonstatespec.part <- t(apply(bincommonmatrix, 1, function(x) { c('0.','1.')[x+1] }))
			colnames(commonstatespec.part) <- intersect.tracks
			commonstatespec.part.reps <- matrix(NA, ncol=length(tracks2use), nrow=nrow(bincommonmatrix))
			for (track in intersect.tracks) {
				index <- which(track==tracks2use)
				commonstatespec.part.reps[,index] <- rep(commonstatespec.part[,as.character(track)], length(index))
			}
			commonstatespec.list[[length(commonstatespec.list)+1]] <- t(apply(commonstatespec.part.reps, 1, function(x) { paste0(x, tracks2use) }))
		}
		commonstatespecs <- do.call(cbind, commonstatespec.list)

		## Go through all commonstate specifications
		binstates.list <- list()
		for (irow in 1:nrow(commonstatespecs)) {
			commonstatespec <- commonstatespecs[irow,]
			commongroups <- levels(factor(commonstatespec))
			binstates.irow <- binstates
			for (commongroup in commongroups) {
				track.index <- which(commonstatespec==commongroup)
				if (grepl('^0\\.', commongroup)) {
					mask <- !apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('|', x) })
				} else if (grepl('^1\\.', commongroup)) {
					mask <- apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('&', x) })
				}
				binstates.irow <- binstates.irow[mask,]
				if (class(binstates.irow)!='matrix') {
					binstates.irow <- matrix(binstates.irow, ncol=length(binstates.irow))
					colnames(binstates.irow) <- tracknames
				}
			}
			binstates.list[[irow]] <- binstates.irow
		}
		binstates.common <- do.call(rbind, binstates.list)
	}
		
	## Merge common and differential states
	if (differential.states | common.states) {
		binstates <- rbind(binstates.diff, binstates.common)
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
