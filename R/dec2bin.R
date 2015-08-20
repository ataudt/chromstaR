#' Conversion of decimal and binary states
#'
#' Convert combinatorial states in decimal representation to combinatorial states in binary representation and vice versa.
#'
#' \pkg{\link{chromstaR}} uses decimal numbers to represent combinatorial states of peaks. These functions serve as a convenient way to get from the efficient decimal representation to a more human-readable binary representation.
#' @return A vector of integers for \code{bin2dec} and a matrix of logicals with one state per row for \code{dec2bin}.
#' @examples
#'## Load example multivariate Hidden Markov Model
#'data(example.multi.HMM)
#'## Get the (decimal) combinatorial states
#'states.decimal <- example.multi.HMM$bins$state
#'## Transform them to binary representation
#'states.binary <- dec2bin(states.decimal,ndigits=length(example.multi.HMM$IDs))
#'## And back to decimal
#'new.states.decimal <- bin2dec(states.binary)
#' @name conversion
NULL
#'
#' @describeIn conversion Decimal to binary conversion.
#' @param dec A vector with whole numbers.
#' @param colnames The column names for the returned matrix. If specified, \code{ndigits} will be the length of \code{colnames}.
#' @param ndigits The number of digits that the binary representation should have. If unspecified, the shortest possible representation will be chosen.
#' @author Aaron Taudt
#' @export
dec2bin = function(dec, colnames=NULL, ndigits=NULL) {

  # Check user input
  dec <- as.numeric(as.character(dec))
	maxdec = max(dec)
	if (!is.null(colnames)) {
		ndigits <- length(colnames)
	} else {
		if (!is.null(ndigits)) {
			if (check.positive.integer(ndigits)!=0) stop("argument 'ndigits' expects a positive integer")
		}
	}
  
  mods <- list()
  mod <- 1
  dectemp <- dec
  while (any(dectemp!=0)) {
    mod <- dectemp %% 2
    mods[[length(mods)+1]] <- mod
    dectemp <- dectemp %/% 2
  }
  
  if (is.null(ndigits)) {
    ndigits <- max(length(mods),1)
  }
  binary.states <- matrix(0, nrow=length(dec), ncol=ndigits)
  if (length(mods)>0) {
    binary.states[,1:ndigits] <- as.matrix(as.data.frame(mods))[,1:ndigits]
  }
  binary.states <- binary.states[,ncol(binary.states):1]
  if (class(binary.states)!='matrix') {
    binary.states <- matrix(binary.states, nrow=length(dec))
  }

	colnames(binary.states) <- colnames
	rownames(binary.states) <- dec
	return(binary.states)

}

#' @describeIn conversion Binary to decimal conversion.
#' @param bin A matrix with only 0 and 1 (or TRUE and FALSE) as entries. One combinatorial state per row.
#' @export
bin2dec = function(bin) {
	if (!is.matrix(bin)) {
		bin <- matrix(bin, nrow=1)
		warning("Argument 'bin' is not a matrix. Interpreting 'bin' as matrix with one row.")
	}
	dec = rep(0,nrow(bin))
	for (i1 in 1:ncol(bin)) {
		dec = dec + 2^(ncol(bin)-i1) * bin[,i1]
	}
	return(dec)
}


