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


