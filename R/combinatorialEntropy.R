#' Combinatorial entropy
#'
#' Calculate the combinatorial entropy along the genome for a set of \code{\link{chromstaR_multivariateHMM}}s.
#'
#' The combinatorial entropy is calculated as ...
#'
#' @author Aaron Taudt
#' @param multi.hmm.list A list of \code{\link{chromstaR_multivariateHMM}} objects or a vector of files that contain such objects.
#' @param direction.of.adding One of \code{c('between.hmms','between.samples')}.
#' @import BiocGenerics
#' @export
combinatorialEntropy <- function(multi.hmm.list, direction.of.adding='between.hmms') {

	if (direction.of.adding == 'between.hmms') {

		message("Calculating entropy ...", appendLF=F); ptm <- proc.time()
		segments <- GRangesList()
		for (i1 in 1:length(multi.hmm.list)) {
			hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[i1]])[[1]] )
			num.bins <- length(hmm$bins)
			num.samples <- length(hmm$IDs)
			binstates <- dec2bin(hmm$segments$state)
			entropy <- apply(binstates, 1, function(x) { log(choose(num.samples, length(which(x)))) })
			isegments <- hmm$segments
			isegments$entropy <- entropy
			segments[[i1]] <- isegments
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		message("Making consensus template ...", appendLF=F); ptm <- proc.time()
		consensus <- GenomicRanges::disjoin(BiocGenerics::unlist(segments))
		conentropy <- matrix(NA, ncol=length(segments), nrow=length(consensus))
		for (i1 in 1:length(segments)) {
			segment <- segments[[i1]]
			mind <- as.matrix(findOverlaps(consensus, segment, select='first'))
			conentropy[,i1] <- segment$entropy[mind]
		}
		mcols(consensus)$entropy <- conentropy
		consensus$sum.entropy <- apply(conentropy, 1, sum)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		# Remove NAs that occur when num.bins differ between HMMs
		consensus <- consensus[!is.na(consensus$sum.entropy)]

		return(consensus)
	}

}


