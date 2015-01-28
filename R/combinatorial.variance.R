combinatorial.variance <- function(multi.hmm, groups, file='combinatorial_variance', window.size.bp=1e5, quantile=0.9) {

	## Transform to GRanges
	if (class(multi.hmm) == class.multivariate.hmm) {
		multi.gr <- multi.hmm$bins
	} else if (class(multi.hmm) == 'GRanges') {
		multi.gr <- multi.hmm
	} else {
		stop("argument 'multi.hmm' expects multivariate HMM object")
	}
	binsize <- width(multi.gr[1])

	## Convert states to binary and split into groups
	binstates <- dec2bin(multi.gr$state, ndigits=ncol(multi.gr$reads))
	binstates.1 <- binstates[,as.logical(groups)]
	binstates.2 <- binstates[,!as.logical(groups)]

	## Variance measure
	diff <- rowSums(abs(binstates.1-binstates.2))
	multi.gr$diff <- diff
	k <- round( window.size.bp / binsize )
	if (k %% 2 == 0) k <- k+1
	var.mean <- as.numeric(runmean(Rle(diff), k, endrule='constant'))
	multi.gr$var.mean <- var.mean

	## Select upper quantile
	cutoff <- quantile(multi.gr$var.mean, quantile)
	mask <- multi.gr$var.mean >= cutoff
	gr <- multi.gr[mask]
	gr <- reduce(gr)

	## Export regions with variance in upper quantile to BED
	# Variables
	filename <- paste0(file,".bed.gz")
	filegz <- gzfile(filename, 'w')
	# Write to file
	message('writing to file',filename)
	cat("", file=filegz)
	numsegments <- length(gr)
	df <- as.data.frame(gr)[,1:3]
	# Adjust coordinates to BED format
	df$start <- df$start - 1
	cat(paste0('track name="combinatorial variance" description="combinatorial variance > ',round(cutoff,3),'" visibility=1 itemRgb=On priority=100\n'), file=filegz, append=TRUE)
	write.table(format(df, scientific=FALSE), file=filegz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	close(filegz)

	## Return modified GRanges
	return(multi.gr)

}
	
