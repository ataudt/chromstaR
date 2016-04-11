#' Normalize read counts
#'
#' Normalize read counts to a given read depth. Reads counts are randomly removed from the input to match the specified read depth.
#'
#' @author Aaron Taudt
#' @param binned.data A \code{\link{GRanges}} object with meta data column 'reads' that contains the read count.
#' @param sample.reads The number of reads that will be retained.
#' @return A \code{\link{GRanges}} object with downsampled read counts.
#' @importFrom stats rbinom
subsample <- function(binned.data, sample.reads) {

	total.reads <- sum(binned.data$counts)
	if (sample.reads >= total.reads) {
		warning("Not resampling reads because sample.reads is bigger than the actual number of reads in the sample.")
		return(binned.data)
	}

	p <- sample.reads / total.reads
	binned.data$counts <- stats::rbinom(binned.data$counts, binned.data$counts, p)

	return(binned.data)
}
	
