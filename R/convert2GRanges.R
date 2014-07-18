binned2GRanges <- function(binned.data, chrom.length.file=NULL, offset=0) {

	library(GenomicRanges)
	gr <- GenomicRanges::GRanges(
			seqnames = Rle(binned.data$chrom),
			ranges = IRanges(start=binned.data$start+offset, end=binned.data$end+offset),
			strand = Rle(strand("*"), nrow(binned.data)),
			reads = binned.data$reads
			)
	if (!is.null(chrom.length.file)) {
		# File with chromosome lengths (1-based)
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		seqlengths(gr) <- as.integer(chrom.lengths[names(seqlengths(gr))])
	}		
	return(gr)

}

hmm2GRanges <- function(hmm, reduce=TRUE) {

# 	library(GenomicRanges)
	### Check user input ###
	if (check.multivariate.model(hmm)!=0 & check.univariate.model(hmm)!=0) stop("argument 'hmm' expects a univariate or multivariate hmm object (type ?uni.hmm or ?multi.hmm for help)")
	if (check.logical(reduce)!=0) stop("argument 'reduce' expects TRUE or FALSE")

	### Create GRanges ###
	# Transfer coordinates
	gr <- GenomicRanges::GRanges(
			seqnames = Rle(hmm$coordinates$chrom),
			ranges = IRanges(start=hmm$coordinates$start, end=hmm$coordinates$end),
			strand = Rle(strand("*"), nrow(hmm$coordinates))
			)
	seqlengths(gr) <- hmm$seqlengths[names(seqlengths(gr))]
	# Reorder seqlevels
	gr <- GenomicRanges::keepSeqlevels(gr, names(hmm$seqlengths))

	if (reduce) {
		# Reduce state by state
		red.gr.list <- GenomicRanges::GRangesList()
		ustates <- unique(hmm$states)
		levels <- levels(hmm$states)
		for (state in ustates) {
			red.gr <- GenomicRanges::reduce(gr[hmm$states==state])
			mcols(red.gr)$states <- rep(factor(state, levels=levels),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		# Merge and sort
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		remove(red.gr.list)
		return(red.gr)
	} else {
		mcols(gr)$reads <- hmm$reads
		mcols(gr)$posteriors <- hmm$posteriors
		mcols(gr)$states <- hmm$states
		return(gr)
	}

}

