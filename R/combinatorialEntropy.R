#' Combinatorial entropy
#'
#' Calculate the combinatorial entropy along the genome for a set of \code{\link{chromstaR_multivariateHMM}}s.
#'
#' The combinatorial entropy describes the variation between samples. For each genomic position, it is defined as
#' \deqn{ log( choose(n,x) ) }
#' where \eqn{n} is the number of samples and \eqn{x} the number of samples that are in state 'modified'. Calculation of this entropy measure only makes sense for samples with the same modification (e.g. 7 samples of H3K36me3 in different tissues). For multiple modifications, the entropy can be calculated for each modification separately and then added.
#'
#' @author Aaron Taudt
#' @param multi.hmm.list A list of \code{\link{chromstaR_multivariateHMM}} objects or a vector of files that contain such objects.
#' @param window.size.bp Window size in base-pairs that will be used to average the results.
#' @param direction.of.adding One of \code{c('between.hmms','inside.hmms')}. Choose 'between.hmms' if your input objects contain combinatorial states from a differential analysis. Choose 'inside.hmms' if your input objects contain combinatorial states from different modifications.
#' @import BiocGenerics
#' @export
combinatorialEntropy <- function(multi.hmm.list, window.size.bp=NULL, direction.of.adding='between.hmms') {

	message("Calculating entropy ...", appendLF=F); ptm <- proc.time()
	if (direction.of.adding == 'between.hmms') {
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
	consensus$entropy <- apply(conentropy, 1, sum)
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	# Remove NAs that occur when num.bins differ between HMMs
	gr <- consensus[!is.na(consensus$entropy)]

	if (!is.null(window.size.bp)) {
		message("Averaging over window size ",window.size.bp,"bp ...", appendLF=F); ptm <- proc.time()
		tg <- unlist(tileGenome(seqlengths(gr), tilewidth=window.size.bp))
		mind <- findOverlaps(gr, tg)
		gr.extended <- gr[queryHits(mind)]
		rlemind <- rle(subjectHits(mind))
		index.last <- cumsum(rlemind$lengths)
		index.first <- c(1,index.last[-length(index.last)]+1)
		start(gr.extended)[index.first] <- start(tg)
		end(gr.extended)[index.last] <- end(tg)
		gr.extended$weighted.entropy <- width(gr.extended)*gr.extended$entropy
		gr.extended$index <- subjectHits(mind)
		df <- as.data.frame(gr.extended)[,c('seqnames','start','end','width','index','weighted.entropy')]
		df <- suppressMessages( collapseBins(df, column2collapseBy='index', columns2sumUp=c('width','weighted.entropy')) )
		gr <- GRanges(seqnames=df$seqnames, ranges=IRanges(start=df$start, end=df$end), entropy=df$sum.weighted.entropy/df$sum.width)
		seqlengths(gr) <- seqlengths(gr.extended)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	}

	return(gr)

}


