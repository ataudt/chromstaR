#' Summary information about called peaks
#'
#' Get information about the peak calls in a \code{\link{multiHMM}} object.
#'
#' @author Aaron Taudt
#' @param multi.hmm A \code{\link{multiHMM}} object.
#' @return A \code{list()} object with various entries.
#' @examples
#'## Load example multivariate Hidden Markov Model
#'data(example.multi.HMM)
#'## Summary information about peaks
#'summarizePeaks(example.multi.HMM)
#' @export
summarizePeaks <- function(multi.hmm) {

	## Get summary statistics about called peaks
	stats <- list()
	stats$total <- list()
	stats$total$ID <- vector()
	stats$total$num.peaks <- vector()
	stats$total$frac.in.peaks <- vector()
	stats$total$mean.len.peaks <- vector()
	stats$total$var.len.peaks <- vector()
	binstates <- dec2bin(multi.hmm$segments$state)
	colnames(binstates) <- multi.hmm$IDs
	segments <- multi.hmm$segments
	for (sample in multi.hmm$IDs) {
		stats$total$ID[sample] <- sample
		stats$total$num.peaks[sample] <- sum(binstates[,sample])
		stats$total$frac.in.peaks[sample] <- sum(as.numeric(width(segments)[binstates[,sample]])) / sum(as.numeric(seqlengths(segments)))
		stats$total$mean.len.peaks[sample] <- mean(width(segments)[binstates[,sample]])
		stats$total$var.len.peaks[sample] <- var(width(segments)[binstates[,sample]])
	}
	stats$total <- as.data.frame(stats$total)
	rownames(stats$total) <- NULL

	## Get summary statistics about called peaks per chromosome
	for (chrom in seqlevels(multi.hmm$bins)) {
		mask <- as.logical(seqnames(multi.hmm$segments)==chrom)
		chrstats <- list()
		chrstats$ID <- vector()
		chrstats$num.peaks <- vector()
		chrstats$frac.in.peaks <- vector()
		chrstats$mean.len.peaks <- vector()
		chrstats$var.len.peaks <- vector()
		chrsegments <- multi.hmm$segments[mask]
		for (sample in multi.hmm$IDs) {
			chrstats$ID[sample] <- sample
			chrstats$num.peaks[sample] <- sum(binstates[mask,sample])
			chrstats$frac.in.peaks[sample] <- sum(as.numeric(width(chrsegments)[binstates[mask,sample]])) / seqlengths(chrsegments)[names(seqlengths(chrsegments))==chrom]
			chrstats$mean.len.peaks[sample] <- mean(width(chrsegments)[binstates[mask,sample]])
			chrstats$var.len.peaks[sample] <- var(width(chrsegments)[binstates[mask,sample]])
		}
		chrstats <- as.data.frame(chrstats)
		rownames(chrstats) <- NULL
		stats[[chrom]] <- chrstats
	}

	## Return
	return(stats)

}

#' Summary information about differential peak calls
#'
#' Get information about the differential peak calls in a \code{\link{multiHMM}} object.
#'
#' @author Aaron Taudt
#' @param multi.hmm A \code{\link{multiHMM}} object.
#' @return A \code{list()} object with various entries.
#' @examples
#'## Load example multivariate Hidden Markov Model
#'data(example.multi.HMM)
#'## Summary information about differentially modified peaks
#'summarizeDiffPeaks(example.multi.HMM)
#' @export
summarizeDiffPeaks <- function(multi.hmm) {

	## Get differential states
	states <- levels(multi.hmm$bins$state)
	binstates <- dec2bin(states)
	diffstates <- states[xor(apply(binstates, 1, function(x) { Reduce('|', x) }), apply(binstates, 1, function(x) { Reduce('&', x) }))]

	diffmask.seg <- multi.hmm$segments$state %in% diffstates
	diffmask.bin <- multi.hmm$bins$state %in% diffstates

	stats <- list()
	stats$num.diff.seg <- vector()
	stats$frac.diff.seg <- vector()
	stats$mean.len.diff.seg <- vector()
	stats$var.len.diff.seg <- vector()

	## Get summary statistics about differentially modified states
	stats$num.diff.seg['total'] <- length(multi.hmm$segments[diffmask.seg])
	stats$frac.diff.seg['total'] <- sum(as.numeric(width(multi.hmm$segments)[diffmask.seg])) / sum(as.numeric(seqlengths(multi.hmm$segments)))
	stats$mean.len.diff.seg['total'] <- mean(width(multi.hmm$segments)[diffmask.seg])
	stats$var.len.diff.seg['total'] <- var(width(multi.hmm$segments)[diffmask.seg])

	## Get summary statistics about differentially modified states per chromosome
	for (chrom in seqlevels(multi.hmm$bins)) {
		stats$num.diff.seg[chrom] <- length(multi.hmm$segments[diffmask.seg & seqnames(multi.hmm$segments)==chrom])
		stats$frac.diff.seg[chrom] <- sum(as.numeric(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])) / sum(as.numeric(seqlengths(multi.hmm$segments)))
		stats$mean.len.diff.seg[chrom] <- mean(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])
		stats$var.len.diff.seg[chrom] <- var(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])
	}
	stats <- as.data.frame(stats)

	## Return
	return(stats)

}

