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
