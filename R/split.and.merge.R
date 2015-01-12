split.per.chrom <- function(uni.hmm, filename=NULL) {

	## Check user input
	if (check.univariate.model(uni.hmm)!=0) {
		uni.hmm <- get(load(uni.hmm))
		if (check.univariate.model(uni.hmm)!=0) stop("argument 'uni.hmm' expects a univariate HMM or a file that contains a univariate HMM")
	}

	## Construct new HMM per chromosome
	split.hmm <- uni.hmm
	split.hmm[c('bins','segments')] <- NULL
	split.bins <- split(uni.hmm$bins, seqnames(uni.hmm$bins))
	split.segments <- split(uni.hmm$segments, seqnames(uni.hmm$segments))

	hmm.list <- list()
	for (chrom in seqlevels(uni.hmm$bins)) {
		chr.hmm <- uni.hmm
		chr.hmm$bins <- split.bins[[chrom]]
		chr.hmm$segments <- split.segments[[chrom]]
		hmm.list[[chrom]] <- chr.hmm
		if (!is.null(filename)) {
			save(chr.hmm, file=paste0(filename,'_chr_',chrom,'.RData'))
		}
	}
	if (is.null(filename)) {
		return(hmm.list)
	}

}

merge.chroms <- function(multi.hmm.list, filename=NULL) {

	## Check user input
	multi.hmm.list <- loadMultiHmmsFromFiles(multi.hmm.list)
		
	## Construct new HMM
	bins <- GRangesList()
	segments <- GRangesList()
	for (i1 in 1:length(multi.hmm.list)) {
		hmm <- multi.hmm.list[[i1]]
		bins[[i1]] <- hmm$bins
		segments[[i1]] <- hmm$segments
	}
	bins <- unlist(bins)
	names(bins) <- NULL
	segments <- unlist(segments)
	names(segments) <- NULL
	multi.hmm <- hmm
	multi.hmm$bins <- bins
	multi.hmm$segments <- segments
	multi.hmm$weights <- table(multi.hmm$bins$state) / length(multi.hmm$bins)

	if (is.null(filename)) {
		return(multi.hmm)
	} else {
		save(multi.hmm, file=filename)
	}

}
