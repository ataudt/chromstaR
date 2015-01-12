summarize <- function(multi.hmm) {

	## Get differential states
	states <- levels(multi.hmm$bins$state)
	binstates <- dec2bin(states)
	diffstates <- states[xor(apply(binstates, 1, function(x) { Reduce('|', x) }), apply(binstates, 1, function(x) { Reduce('&', x) }))]

	diffmask.seg <- multi.hmm$segments$state %in% diffstates
	diffmask.bin <- multi.hmm$bins$state %in% diffstates

	## Get summary statistics about differentially modified states
	stats <- list()
	stats$num.diff.seg <- length(multi.hmm$segments[diffmask.seg])
	stats$frac.diff.seg <- sum(as.numeric(width(multi.hmm$segments)[diffmask.seg])) / sum(as.numeric(seqlengths(multi.hmm$segments)))
	stats$mean.len.diff.seg <- mean(width(multi.hmm$segments)[diffmask.seg])
	stats$var.len.diff.seg <- var(width(multi.hmm$segments)[diffmask.seg])

	## Get summary statistics about differentially modified states per chromosome
	stats$diff.per.chrom$num.diff.seg <- vector()
	stats$diff.per.chrom$frac.diff.seg <- vector()
	stats$diff.per.chrom$mean.len.diff.seg <- vector()
	stats$diff.per.chrom$var.len.diff.seg <- vector()
	for (chrom in seqlevels(multi.hmm$bins)) {
		stats$diff.per.chrom$num.diff.seg[chrom] <- length(multi.hmm$segments[diffmask.seg & seqnames(multi.hmm$segments)==chrom])
		stats$diff.per.chrom$frac.diff.seg[chrom] <- sum(as.numeric(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])) / sum(as.numeric(seqlengths(multi.hmm$segments)))
		stats$diff.per.chrom$mean.len.diff.seg[chrom] <- mean(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])
		stats$diff.per.chrom$var.len.diff.seg[chrom] <- var(width(multi.hmm$segments)[diffmask.seg & as.logical(seqnames(multi.hmm$segments)==chrom)])
	}
	stats$diff.per.chrom <- as.data.frame(stats$diff.per.chrom)

	## Get summary statistics per track
	stats$per.track <- list()
	stats$per.track$ID <- vector()
	stats$per.track$num.peaks <- vector()
	stats$per.track$frac.in.peaks <- vector()
	stats$per.track$mean.len.peaks <- vector()
	stats$per.track$var.len.peaks <- vector()
	binstates <- dec2bin(multi.hmm$segments$state)
	colnames(binstates) <- multi.hmm$IDs.univariate
	for (track in multi.hmm$IDs.univariate) {
		stats$per.track$ID[track] <- track
		stats$per.track$num.peaks[track] <- sum(binstates[,track])
		stats$per.track$frac.in.peaks[track] <- sum(as.numeric(width(multi.hmm$segments)[binstates[,track]])) / sum(as.numeric(seqlengths(multi.hmm$segments)))
		stats$per.track$mean.len.peaks[track] <- mean(width(multi.hmm$segments)[binstates[,track]])
		stats$per.track$var.len.peaks[track] <- var(width(multi.hmm$segments)[binstates[,track]])
	}
	stats$per.track <- as.data.frame(stats$per.track)
	rownames(stats$per.track) <- NULL

	## Return
	return(stats)

}
