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
			mcols(red.gr)$state <- rep(factor(state, levels=levels),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		# Merge and sort
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		remove(red.gr.list)
		return(red.gr)
	} else {
		mcols(gr)$reads <- hmm$reads
		mcols(gr)$posteriors <- hmm$posteriors
		mcols(gr)$state <- hmm$states
		return(gr)
	}

}

bed2GRanges <- function(bedfile, chrom.length.file, skip=1, binsize=NULL) {

	# File with chromosome lengths (1-based)
	chrom.lengths.df <- read.table(chrom.length.file)
	chrom.lengths <- chrom.lengths.df[,2]
	names(chrom.lengths) <- chrom.lengths.df[,1]
	# File with reads, determine classes first for faster import (0-based)
	tab5rows <- read.table(bedfile, nrows=5, skip=skip)
	classes.in.bed <- sapply(tab5rows, class)
	classes <- rep("NULL",length(classes.in.bed))
	classes[1:4] <- classes.in.bed[1:4]
	data <- read.table(bedfile, colClasses=classes, skip=skip)
	# Convert to GRanges object
	data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]),
																	ranges=IRanges(start=data[,2]+1, end=data[,3]+1),	# +1 to match coordinate systems
																	strand=Rle(strand("*"), nrow(data)),
																	state=data[,4])
	seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	chroms.in.data <- seqlevels(data)

# Also possible but slower!
# 	## Create binned GRanges object
# 	binned <- GenomicRanges::GRanges()
# 	for (chrom in mixedsort(chroms.in.data)) {
# 		# Check last incomplete bin
# 		incomplete.bin <- chrom.lengths[chrom] %% binsize > 0
# 		if (incomplete.bin) {
# 			numbins <- floor(chrom.lengths[chrom]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
# 		} else {
# 			numbins <- chrom.lengths[chrom]/binsize
# 		}
# 		ibinned <- GenomicRanges::GRanges(seqnames=Rle(rep(chrom, numbins)),
# 																				ranges=IRanges(start=seq(from=1, by=binsize, length=numbins), end=seq(from=binsize, by=binsize, length=numbins)))
# 		seqlengths(ibinned) <- chrom.lengths[chrom]
# 		suppressWarnings( binned <- c(binned, ibinned) )
# 	}
# 	mcols(binned)$state <- NA
# 
# 	## Find overlaps and assign states
# 	data.state <- split(data, mcols(data)$state)
# 	mind <- as.matrix(findOverlaps(binned, data.state))	# subjectHits is the name of the list entry, which is the state
# 	mcols(binned)$state[mind[,'queryHits']] <- mind[,'subjectHits']
# 	mcols(data)$state <- factor(mcols(data)$state)
# 
# 	return(binned)

	## Inflate every range with bins
	if (!is.null(binsize)) {
		grl <- split(data, seqnames(data))
		inflated.data <- GRangesList()
		for (i1 in 1:length(grl)) {
			rgr <- ranges(grl[[i1]])
			widths <- width(rgr)
			numbins <- widths %/% binsize
			starts <- start(rgr)
			ends <- end(rgr)
			chroms <- seqnames(grl[[i1]])
			states <- mcols(grl[[i1]])$state

			# Create inflated vectors
			rle <- rle(1)
			rle$lengths <- numbins
			rle$values <- as.character(chroms)
			infchroms <- inverse.rle(rle)
			rle$values <- states
			infstates <- inverse.rle(rle)
			infstarts <- seq(starts[1], ends[length(ends)]-1, by=binsize)
			infends <- seq(starts[1]-1+binsize, ends[length(ends)]-2+binsize, by=binsize)

			inflated.chrom <- GenomicRanges::GRanges(seqnames=Rle(infchroms),
																								ranges=IRanges(start=infstarts, end=infends),
																								strand=Rle(strand('*'), sum(numbins)),
																								state=infstates)
			suppressWarnings( inflated.data[[i1]] <- inflated.chrom )
		}
		inflated.data <- unlist(inflated.data)
		seqlengths(inflated.data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		data <- inflated.data
	}
	mcols(data)$state <- factor(mcols(data)$state)

	return(data)

}


