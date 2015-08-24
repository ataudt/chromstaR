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
#' @param direction.of.adding Calculation of the combinatorial entropy only makes sense between different modifications. Choose 'between.hmms' if you have one modification per HMM, which is the case in a differential analysis. Choose 'inside.hmms' if you have multiple modifications per HMM, which is the case in a combinatorial analysis. In this case it is up to the user to make sure that the \code{IDs} between HMMs are identical.
#' @import BiocGenerics
#' @export
combinatorialEntropy <- function(multi.hmm.list, window.size.bp=NULL, direction.of.adding='between.hmms') {

	if (direction.of.adding == 'between.hmms') {
		#===================================
		### Direction of adding between HMMs
		#===================================
		message("Calculating entropy ...", appendLF=F); ptm <- proc.time()
		segments <- GRangesList()
		for (i1 in 1:length(multi.hmm.list)) {
			hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[i1]])[[1]] )
			num.samples <- length(hmm$IDs)
			binstates <- dec2bin(hmm$segments$state)
			num.1s <- rowSums(binstates)
			entropy <- log(choose(num.samples, num.1s))
			isegments <- hmm$segments
			mcols(isegments) <- NULL
			isegments$entropy <- entropy
			segments[[i1]] <- isegments
		}
		remove(hmm)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		message("Making consensus template ...", appendLF=F); ptm <- proc.time()
		consensus <- GenomicRanges::disjoin(BiocGenerics::unlist(segments))
		conentropy <- matrix(NA, ncol=length(segments), nrow=length(consensus))
		for (i1 in 1:length(segments)) {
			segment <- segments[[i1]]
			mind <- as.matrix(findOverlaps(consensus, segment, select='first'))
			conentropy[,i1] <- segment$entropy[mind]
		}
		consensus$entropy <- rowSums(conentropy)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	} else if (direction.of.adding == 'inside.hmms') {
		#==================================
		### Direction of adding inside HMMs
		#==================================
		message("Making consensus templates ...", appendLF=F); ptm <- proc.time()
		segments <- GRangesList()
		binstates <- list()
		for (i1 in 1:length(multi.hmm.list)) {
			hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[i1]])[[1]] )
			segments[[i1]] <- hmm$segments[,'state']
			binstates[[i1]] <- dec2bin(hmm$segments$state)
			colnames(binstates[[i1]]) <- hmm$IDs
		}
		remove(hmm)
		consensus <- GenomicRanges::disjoin(BiocGenerics::unlist(segments))
		IDs <- sort(unique(unlist(lapply(binstates,colnames))))
		conbinstates <- array(NA, dim=c(length(consensus), length(segments), length(IDs)), dimnames=list(region=1:length(consensus), model=1:length(segments), ID=IDs))
		for (i1 in 1:length(segments)) {
			segment <- segments[[i1]]
			mind <- as.matrix(findOverlaps(consensus, segment, select='first'))
			for (ID in IDs) {
				colindex <- which(colnames(binstates[[i1]])==ID)
				if (length(colindex)!=0) {
					conbinstates[,i1,ID] <- binstates[[i1]][mind,colindex]
				}
			}
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		message("Calculating entropy ...", appendLF=F); ptm <- proc.time()
		entropy <- matrix(NA, nrow=length(consensus), ncol=length(IDs))
		colnames(entropy) <- IDs
		for (i1 in 1:length(IDs)) {
			num.samples <- length(which(!is.na(conbinstates[1,,i1])))
			num.1s <- rowSums(conbinstates[,,i1], na.rm=T)
			entropy[,i1] <- log(choose(num.samples, num.1s))
		}
		consensus$entropy <- rowSums(entropy)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	}

	#===========================
	### Average over window size
	#===========================
	# Remove NAs that occur when num.bins differ between HMMs
	gr <- consensus[!is.na(consensus$entropy)]
	if (!is.null(window.size.bp)) {
		message("Averaging over window size ",window.size.bp,"bp ...", appendLF=F); ptm <- proc.time()
		tg <- unlist(tileGenome(seqlengths(gr), tilewidth=window.size.bp))
		tg.cut <- subsetByOverlaps(tg, gr)
		mind <- findOverlaps(gr, tg.cut)
		gr.extended <- gr[queryHits(mind)]
		rlemind <- rle(subjectHits(mind))
		index.last <- cumsum(rlemind$lengths)
		index.first <- c(1,index.last[-length(index.last)]+1)
		start(gr.extended)[index.first] <- start(tg.cut)
		end(gr.extended)[index.last] <- end(tg.cut)
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


