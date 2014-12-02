dec2bin = function(dec, ndigits=NULL, names=NULL) {

	# Convert factor to integer
	dec <- as.integer(as.character(dec))
	# Check user input
	maxdec = max(dec)
	if (is.null(ndigits)) {
		if (is.null(names)) {
			ndigits = max(which(as.logical(intToBits(maxdec))))
		} else {
			ndigits = length(names)
		}
	}
	if (!is.null(names)) {
		if (maxdec >= 2^length(names)) {
			stop("Cannot transform state ",max(dec)," into name representation. Cause: Not enough names or state to high.")
		}
	}

	binary_states = matrix(as.logical(intToBits(dec)), nrow=length(dec), byrow=TRUE)
	binary_states = binary_states[ ,ndigits:1]

	if (!is.null(names)) {
		namemat = matrix(rep(names, length(dec)), nrow=length(dec), byrow=TRUE)
		namemat[!binary_states] = ""
		dfnames = data.frame(
			comb.state = dec,
			namemat
		)
		names(dfnames)[-1] = names
		return(dfnames)
	} else {
		return(binary_states)
	}

}

bin2dec = function(bin) {
	dec = rep(0,nrow(bin))
	for (i1 in 1:ncol(bin)) {
		dec = dec + 2^(ncol(bin)-i1) * bin[,i1]
	}
	return(dec)
}

state.brewer <- function(statespec) {

	# List of possible binary states
	ntracks <- length(statespec)
	all.binstates <- dec2bin(0:(2^ntracks-1))

	# Generate select mask
	mask <- rep(TRUE, nrow(all.binstates))
	for (itrack in 1:length(statespec)) {
		bit <- statespec[itrack]
		if (bit == 1 | bit == 0) {
			mask <- mask & all.binstates[,itrack]==bit
		}
	}

	# Select specified states
	binstates <- all.binstates[mask,]
	combstates <- bin2dec(binstates)

	return(combstates)

}
