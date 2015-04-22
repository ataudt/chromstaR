#' Enrichment of (combinatorial) states for genomic annotations
#'
#' The function calculates the enrichment of a genomic feature with peaks or combinatorial states. Input is a \code{\link{chromstaR_multivariateHMM}} object (containing the peak calls and combinatorial states) and a \code{\link{GRanges}} object containing the annotation of interest (e.g. transcription start sites or genes).
#'
#' @author Aaron Taudt
#' @param hmm A \code{\link{chromstaR_multivariateHMM}} object.
#' @param annotation A \code{\link{GRanges}} object with the annotation of interest.
#' @param bp.around.annotation An integer specifying the number of basepairs up- and downstream of the annotation for which the enrichment will be calculated.
#' @param region A combination of \code{c('start','inside','end')} specifying the region of the annotation for which the enrichment will be calculated.
#' @return A \code{list()} containing \code{data.frame()}s for enrichment of combinatorial states, binary states and posterior probabilities at the start, end and inside of the annotation.
#' @import S4Vectors
#' @export
enrichmentCurve <- function(hmm, annotation, bp.around.annotation=10000, region=c('start','inside','end')) {

	## Variables
	binsize <- width(hmm$bins)[1]
	lag <- round(bp.around.annotation/binsize)
	eCurve <- list()
	eCurve$combstates <- list()
	eCurve$binstates <- list()
	eCurve$posteriors <- list()
	eCurve$reads <- list()

	## Get combinatorial and binary states
	combstates <- hmm$bins$state
	binstates <- dec2bin(hmm$bins$state, colnames=hmm$IDs)

	### Calculationg enrichment inside of annotation ###
	if ('inside' %in% region) {
		message("Enrichment inside of annotations")

		# Get bins that overlap annotation
		ind <- findOverlaps(hmm$bins, annotation)
		bins.per.annotation <- table(subjectHits(ind))	# Table is sorted!
		annotation$num.bins.spanning <- rep(0, length(annotation))
		annotation$num.bins.spanning[as.numeric(names(bins.per.annotation))] <- bins.per.annotation
		strand.per.annotation <- S4Vectors::as.factor(strand(annotation[as.numeric(names(bins.per.annotation))]))
		names(strand.per.annotation) <- names(bins.per.annotation)
		
		# States, posteriors and strands per bin-that-overlaps-an-annotation
		anno.binstates <- binstates[queryHits(ind),]
		anno.combstates <- combstates[queryHits(ind)]
		anno.post <- hmm$bins$posteriors[queryHits(ind),]
		colnames(anno.post) <- NULL
		anno.reads <- hmm$bins$reads[queryHits(ind),]
		colnames(anno.reads) <- NULL
		anno.strands <- strand(annotation)[subjectHits(ind)]

		# Relative coordinate of every bin (TSS=0, TTS=1)
		relcoord <- list()
		for (i1 in 1:length(bins.per.annotation)) {
			if (strand.per.annotation[i1]=='+' | strand.per.annotation[i1]=='*') {
				relcoord[[length(relcoord)+1]] <- (1:bins.per.annotation[i1] - 1) / (bins.per.annotation[i1]-1)
			} else if (strand.per.annotation[i1]=='-') {
				relcoord[[length(relcoord)+1]] <- rev( (1:bins.per.annotation[i1] - 1) / (bins.per.annotation[i1]-1) )
			}
		}
		relcoord <- unlist(relcoord)

		# Collect in data.frame
		anno.df <- data.frame(as.data.frame(ind), strand=anno.strands, binstate=anno.binstates, post=anno.post, reads=anno.reads, combstate=anno.combstates)
		# Reorder to add stuff that was calculated from sorted table
		anno.df <- cbind(anno.df[order(anno.df$subjectHits),], relcoord=relcoord)
		# Annotations that only fall into 1 bin need to get both interval 0 (start) and 1 (end)
		temp <- anno.df[is.nan(anno.df$relcoord),]
		anno.df <- anno.df[!is.nan(anno.df$relcoord),]
		anno.df <- rbind(anno.df, temp)
		anno.df$relcoord[is.nan(anno.df$relcoord)] <- 0
		anno.df <- rbind(anno.df, temp)
		anno.df$relcoord[is.nan(anno.df$relcoord)] <- 1

		# Interval
		intervals <- sort(c(-0.1, seq(from=1e-9, to=1, length=11), 1.1))
		intervals <- sort(seq(from=0, to=1, length=101))
		anno.df$interval <- intervals[findInterval(anno.df$relcoord, intervals)]

		# Mean over intervals
		binstates.inside <- matrix(NA, nrow=length(intervals), ncol=ncol(anno.binstates), dimnames=list(interval=intervals, track=hmm$IDs))
		combstates.inside <- matrix(NA, nrow=length(intervals), ncol=length(levels(anno.df$combstate)), dimnames=list(interval=intervals, state=levels(anno.df$combstate)))
		posteriors.inside <- matrix(NA, nrow=length(intervals), ncol=ncol(anno.binstates), dimnames=list(interval=intervals, track=hmm$IDs))
		reads.inside <- matrix(NA, nrow=length(intervals), ncol=ncol(anno.binstates), dimnames=list(interval=intervals, track=hmm$IDs))
		for (interval in intervals) {
			i1 <- which(interval==intervals)
			mask <- anno.df$interval==interval
			binstates.inside[i1,] <- colMeans(anno.df[,grepl('binstate',names(anno.df))][mask,], na.rm=T)
			posteriors.inside[i1,] <- colMeans(anno.df[,grepl('post',names(anno.df))][mask,], na.rm=T)
			reads.inside[i1,] <- colMeans(anno.df[,grepl('reads',names(anno.df))][mask,], na.rm=T)
			combstates.inside[i1,] <- table(anno.df[grepl('combstate',names(anno.df))][mask,])
		}

		rownames(combstates.inside) <- as.numeric(rownames(combstates.inside)) * binsize
		rownames(binstates.inside) <- as.numeric(rownames(binstates.inside)) * binsize
		rownames(posteriors.inside) <- as.numeric(rownames(posteriors.inside)) * binsize
		rownames(reads.inside) <- as.numeric(rownames(reads.inside)) * binsize
		eCurve$combstates$inside <- combstates.inside
		eCurve$binstates$inside <- binstates.inside
		eCurve$posteriors$inside <- posteriors.inside
		eCurve$reads$inside <- reads.inside
	}
	
	### 10000 bp before and after annotation ###
	if ('start' %in% region) {
		message("Enrichment ",bp.around.annotation,"bp before annotations")
		# Get bins that overlap the start of annotation
		index.start.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], hmm$bins, select="first")
		index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="last")
		index.start.plus <- index.start.plus[!is.na(index.start.plus)]
		index.start.minus <- index.start.minus[!is.na(index.start.minus)]
		# Occurrences at every bin position relative to feature
		binstates.start <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		combstates.start <- array(dim=c(length(-lag:lag), length(levels(hmm$bins$state))), dimnames=list(lag=-lag:lag, state=levels(hmm$bins$state)))
		posteriors.start <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		reads.start <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		for (ilag in -lag:lag) {
	# 		message("lag = ",ilag)
			index <- c(index.start.plus+ilag, index.start.minus-ilag)
			index <- index[index>0 & index<=nrow(binstates)]
			binstates.start[as.character(ilag),] <- colMeans(binstates[index,])
			combstates.start[as.character(ilag),] <- table(combstates[index])
			posteriors.start[as.character(ilag),] <- colMeans(hmm$bins$posteriors[index,])
			reads.start[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		}
		rownames(combstates.start) <- as.numeric(rownames(combstates.start)) * binsize
		rownames(binstates.start) <- as.numeric(rownames(binstates.start)) * binsize
		rownames(posteriors.start) <- as.numeric(rownames(posteriors.start)) * binsize
		rownames(reads.start) <- as.numeric(rownames(reads.start)) * binsize
		eCurve$combstates$start <- combstates.start
		eCurve$binstates$start <- binstates.start
		eCurve$posteriors$start <- posteriors.start
		eCurve$reads$start <- reads.start
	}
	if ('end' %in% region) {
		message("Enrichment ",bp.around.annotation,"bp after annotations")
		# Get bins that overlap the end of annotation
		index.end.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], hmm$bins, select="last")
		index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="first")
		index.end.plus <- index.end.plus[!is.na(index.end.plus)]
		index.end.minus <- index.end.minus[!is.na(index.end.minus)]
		# Occurrences at every bin position relative to feature
		binstates.end <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		combstates.end <- array(dim=c(length(-lag:lag), length(levels(hmm$bins$state))), dimnames=list(lag=-lag:lag, state=levels(hmm$bins$state)))
		posteriors.end <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		reads.end <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		for (ilag in -lag:lag) {
	# 		message("lag = ",ilag)
			index <- c(index.end.plus+ilag, index.end.minus-ilag)
			index <- index[index>0 & index<=nrow(binstates)]
			binstates.end[as.character(ilag),] <- colMeans(binstates[index,])
			combstates.end[as.character(ilag),] <- table(combstates[index])
			posteriors.end[as.character(ilag),] <- colMeans(hmm$bins$posteriors[index,])
			reads.end[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		}
		rownames(combstates.end) <- as.numeric(rownames(combstates.end)) * binsize
		rownames(binstates.end) <- as.numeric(rownames(binstates.end)) * binsize
		rownames(posteriors.end) <- as.numeric(rownames(posteriors.end)) * binsize
		rownames(reads.end) <- as.numeric(rownames(reads.end)) * binsize
		eCurve$combstates$end <- combstates.end
		eCurve$binstates$end <- binstates.end
		eCurve$posteriors$end <- posteriors.end
		eCurve$reads$end <- reads.end
	}

	return(eCurve)

}

plot.cross.correlation <- function(hmm, annotation, bp.around.annotation=10000) {

	## Debugging
	library(biomaRt)
	bp.around.annotation=10000
# 	hg19 <- useMart('ENSEMBL_MART_ENSEMBL', host='grch37.ensembl.org', dataset='hsapiens_gene_ensembl')
# 	filters <- listFilters(hg19)
# 	attributes <- listAttributes(hg19)
# 	hg19.genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'), mart=hg19)
# 	names(hg19.genes)[1:5] <- c('ID','chrom','start','end','strand')
# 	genes <- GRanges(seqnames=paste0('chr',hg19.genes$chrom), ranges=IRanges(start=hg19.genes$start, end=hg19.genes$end), strand=hg19.genes$strand, gene_id=hg19.genes$ID)
# 	file <- "multiresults_all_chrom/multivariate_mark_H3K27ac_patient_149_binsize_1000.RData"
# 	hmm <- get(load(file))
# 	annotation <- genes

	rat.ensembl59 <- useMart('ENSEMBL_MART_ENSEMBL', host='aug2010.archive.ensembl.org', dataset='rnorvegicus_gene_ensembl')
	filters <- listFilters(rat.ensembl59)
	attributes <- listAttributes(rat.ensembl59)
	rat.genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'), mart=rat.ensembl59)
	names(rat.genes)[1:5] <- c('ID','chrom','start','end','strand')
	genes <- GRanges(seqnames=paste0('chr',rat.genes$chrom), ranges=IRanges(start=rat.genes$start, end=rat.genes$end), strand=rat.genes$strand, gene_id=rat.genes$ID)
	file <- 'multiresults/euratrans_lv_binsize_1000.RData'
	hmm <- get(load(file))
	annotation <- genes

	## Variables
	binsize <- width(hmm$bins)[1]
	lag <- round(bp.around.annotation/binsize)

	### Calculationg cross-correlation inside of normed annotation ###
	message("calculating cross-correlation inside of annotation")
	# Select only annotations that span more than 10 bins
	annotation.sub <- annotation
	# Not necessary but makes it easier to debug
	annotation.sub <- keepSeqlevels(annotation.sub, seqlevels(hmm$bins))
	seqlengths(annotation.sub) <- seqlengths(hmm$bins)
	annotation.sub <- sort(annotation.sub)

	# Get bins that overlap annotation
	ind <- findOverlaps(hmm$bins, annotation.sub)
	bins.per.annotation <- table(subjectHits(ind))	# Table is sorted!
	annotation.sub$num.bins.spanning <- rep(0, length(annotation.sub))
	annotation.sub$num.bins.spanning[as.numeric(names(bins.per.annotation))] <- bins.per.annotation
	strand.per.annotation <- as.factor(strand(annotation.sub[as.numeric(names(bins.per.annotation))]))
	names(strand.per.annotation) <- names(bins.per.annotation)
	
	# Get binary states
	binstates <- dec2bin(hmm$bins$state, ndigits=ncol(hmm$bins$reads))

	# States, reads and strands per bin-that-overlaps-an-annotation
	ind.binstates <- binstates[queryHits(ind),]
	ind.reads <- hmm$bins$reads[queryHits(ind),]
	colnames(ind.reads) <- NULL
	ind.strands <- strand(annotation.sub)[subjectHits(ind)]

	# Relative coordinate of every bin (TSS=0, TTS=1)
	relcoord <- list()
	for (i1 in 1:length(bins.per.annotation)) {
		if (strand.per.annotation[i1]=='+') {
			relcoord[[length(relcoord)+1]] <- (1:bins.per.annotation[i1] - 1) / (bins.per.annotation[i1]-1)
		} else if (strand.per.annotation[i1]=='-') {
			relcoord[[length(relcoord)+1]] <- rev( (1:bins.per.annotation[i1] - 1) / (bins.per.annotation[i1]-1) )
		}
	}
	ind.relcoord <- unlist(relcoord)

	# Collect in data.frame
	ind.df <- data.frame(as.data.frame(ind), strand=ind.strands, binstate=ind.binstates, read=ind.reads)
	# Reorder to add stuff that was calculated from sorted table
	ind.df <- cbind(ind.df[order(ind.df$subjectHits),], relcoord=ind.relcoord)
	# Annotations that only fall into 1 bin need to get both interval 0 (start) and 1 (end)
	temp <- ind.df[is.nan(ind.df$relcoord),]
	ind.df <- ind.df[!is.nan(ind.df$relcoord),]
	ind.df <- rbind(ind.df, temp)
	ind.df$relcoord[is.nan(ind.df$relcoord)] <- 0
	ind.df <- rbind(ind.df, temp)
	ind.df$relcoord[is.nan(ind.df$relcoord)] <- 1

	# Interval
	intervals <- sort(c(-0.1, seq(from=1e-9, to=1, length=11), 1.1))
	intervals <- sort(seq(from=0, to=1, length=101))
	ind.df$interval <- intervals[findInterval(ind.df$relcoord, intervals)]

	# Mean over intervals
	binstates.inside <- matrix(NA, nrow=length(intervals), ncol=ncol(ind.binstates), dimnames=list(interval=intervals, track=colnames(hmm$bins$reads)))
	reads.inside <- matrix(NA, nrow=length(intervals), ncol=ncol(ind.binstates), dimnames=list(interval=intervals, track=colnames(hmm$bins$reads)))
	for (interval in intervals) {
		i1 <- which(interval==intervals)
		mask <- ind.df$interval==interval
		binstates.inside[i1,] <- colMeans(ind.df[,grepl('binstate',names(ind.df))][mask,], na.rm=T)
		reads.inside[i1,] <- colMeans(ind.df[,grepl('read',names(ind.df))][mask,], na.rm=T)
	}
	
	### 10000 bp before and after annotation ###
	message("calculating cross-correlation before annotation")
	# Get bins that overlap the start of annotation
	index.start.plus <- findOverlaps(annotation[strand(annotation)=='+'], hmm$bins, select="first")
	index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="last")
	index.start.plus <- index.start.plus[!is.na(index.start.plus)]
	index.start.minus <- index.start.minus[!is.na(index.start.minus)]
	# Occurrences at every bin position relative to feature
	binstates.before <- array(dim=c(length(-lag:0), ncol(hmm$bins$reads)), dimnames=list(lag=-lag:0, track=colnames(hmm$bins$reads)))
	reads.before <- array(dim=c(length(-lag:0), ncol(hmm$bins$reads)), dimnames=list(lag=-lag:0, track=colnames(hmm$bins$reads)))
	for (ilag in -lag:0) {
		message("lag = ",ilag,"\r", appendLF=F)
		index <- c(index.start.plus+ilag, index.start.minus-ilag)
		index <- index[index>0]
		binstates.before[as.character(ilag),] <- colMeans(binstates[index,])
		reads.before[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		message("             \r")
	}
	message("calculating cross-correlation after annotation")
	# Get bins that overlap the end of annotation
	index.end.plus <- findOverlaps(annotation[strand(annotation)=='+'], hmm$bins, select="last")
	index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="first")
	index.end.plus <- index.end.plus[!is.na(index.end.plus)]
	index.end.minus <- index.end.minus[!is.na(index.end.minus)]
	# Occurrences at every bin position relative to feature
	binstates.after <- array(dim=c(length(0:lag), ncol(hmm$bins$reads)), dimnames=list(lag=0:lag, track=colnames(hmm$bins$reads)))
	reads.after <- array(dim=c(length(0:lag), ncol(hmm$bins$reads)), dimnames=list(lag=0:lag, track=colnames(hmm$bins$reads)))
	for (ilag in 0:lag) {
		message("lag = ",ilag,"\r", appendLF=F)
		index <- c(index.end.plus+ilag, index.end.minus-ilag)
		index <- index[index>0]
		binstates.after[as.character(ilag),] <- colMeans(binstates[index,])
		reads.after[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		message("             \r", appendLF=F)
	}

	### Combine results ###
	positions <- c(-lag:0 * binsize, intervals, 0:lag * binsize)
	positions.plot <- c(-lag:0 * binsize, (intervals * lag * binsize), lag*binsize + 0:lag * binsize)

	binstates.around <- as.data.frame(rbind(binstates.before, binstates.inside, binstates.after))
	rownames(binstates.around) <- NULL
	binstates.around$positions <- positions
	binstates.around$positions.plot <- positions.plot
	reads.around <- as.data.frame(rbind(reads.before, reads.inside, reads.after))
	reads.around$positions <- positions
	reads.around$positions.plot <- positions.plot
	
	df <- melt(reads.around, varnames=c('pos','track'), id.vars=c('positions','positions.plot'), variable.name='track')
	ggplt.reads <- ggplot(df) + geom_line(aes_string(x='positions.plot', y='value', col='track'))
	df <- melt(binstates.around, varnames=c('pos','track'), id.vars=c('positions','positions.plot'), variable.name='track')
	ggplt.binstates <- ggplot(df) + geom_line(aes_string(x='positions.plot', y='value', col='track'))


	
}

cross.correlation <- function(multi.hmm, annotation.file, grouping=NULL, lag.in.bp=10000) {

	## Import annotation file
	message("importing ",annotation.file," ...", appendLF=F)
	anno.data <- rtracklayer::import(annotation.file)
	message(" done")

	## Convert multi.hmm to GRanges
	if (class(multi.hmm) == class.multivariate.hmm) {
		message("converting multi.hmm to GRanges...", appendLF=F)
		multi.gr <- multi.hmm$bins
		message(" done")
	} else if (class(multi.hmm) == 'GRanges') {
		multi.gr <- multi.hmm
	} else {
		stop("argument 'multi.hmm' expects multivariate HMM object")
	}

	## Calculate variables
	features <- levels(anno.data$type)
	binsize <- width(multi.gr)[1]
	lag <- round(lag.in.bp/binsize)
	if (is.null(grouping)) {
		grouping <- rep(0,ncol(multi.gr$reads))
	}

	## Split the combinatorial states per group
	states.per.group <- list()
	binstates <- dec2bin(multi.gr$state, ndigits=ncol(multi.gr$reads))
	for (group in grouping) {
		num.group.members <- length(which(group==grouping))
		states.per.group[[as.character(group)]] <- factor(bin2dec(binstates[,group==grouping]), levels=0:(2^num.group.members-1))
	}
	groups <- unique(grouping)
	levels.group <- unique(unlist(lapply(states.per.group, levels)))

	### Loop over every feature in the data
	ggplts <- list()
	for (feature in features) {
		message("Calculating cross correlation for feature '",feature,"'")
		anno.data.feature <- anno.data[anno.data$type==feature]
		## Do stuff for each strand separately to take directionality of the features into account
		# Get bins that overlap the start of that feature
		index.start.plus <- findOverlaps(anno.data.feature[strand(anno.data.feature)=='+'], multi.gr, select="first")
		index.start.minus <- findOverlaps(anno.data.feature[strand(anno.data.feature)=='-'], multi.gr, select="last")
		index.start.plus <- index.start.plus[!is.na(index.start.plus)]
		index.start.minus <- index.start.minus[!is.na(index.start.minus)]
		# Get bins that overlap the end of that feature
		index.end.plus <- findOverlaps(anno.data.feature[strand(anno.data.feature)=='+'], multi.gr, select="last")
		index.end.minus <- findOverlaps(anno.data.feature[strand(anno.data.feature)=='-'], multi.gr, select="first")
		index.end.plus <- index.end.plus[!is.na(index.end.plus)]
		index.end.minus <- index.end.minus[!is.na(index.end.minus)]
		# Occurrences at every bin position relative to feature
		tables.start <- array(dim=c(length(-lag:lag), length(levels.group), length(groups)), dimnames=list(lag=-lag:lag, comb.state=levels.group, group=groups))
		tables.end <- tables.start
		tracks.start <- array(dim=c(length(-lag:lag), ncol(multi.gr$reads)), dimnames=list(lag=-lag:lag, track=colnames(multi.gr$reads)))
		tracks.end <- tracks.start
		for (ilag in -lag:lag) {
			message("lag = ",ilag,"\r", appendLF=F)
			# Combinatorial states
			for (group in groups) {
				tables.start[as.character(ilag), ,as.character(group)] <- table(states.per.group[[as.character(group)]][index.start.plus+ilag]) + table(states.per.group[[as.character(group)]][index.start.minus-ilag])
				tables.end[as.character(ilag), ,as.character(group)] <- table(states.per.group[[as.character(group)]][index.end.plus+ilag]) + table(states.per.group[[as.character(group)]][index.end.minus-ilag])
			}
			# Tracks
			tracks.start[as.character(ilag),] <- colSums(dec2bin(multi.gr[index.start.plus+ilag]$state, ndigits=ncol(multi.gr$reads))) + colSums(dec2bin(multi.gr[index.start.minus-ilag]$state, ndigits=ncol(multi.gr$reads)))
			tracks.end[as.character(ilag),] <- colSums(dec2bin(multi.gr[index.end.plus+ilag]$state, ndigits=ncol(multi.gr$reads))) + colSums(dec2bin(multi.gr[index.end.minus-ilag]$state, ndigits=ncol(multi.gr$reads)))
			message("             \r", appendLF=F)
		}

		## Plot cross-correlations
		# Combinatorial states
		df <- melt(tables.start)
		df$lag.in.bp <- df$lag * binsize
		df$comb.state <- factor(df$comb.state, levels=levels(multi.gr$state))
		df$group <- factor(df$group, levels=levels.group)
		ggplts[[as.character(feature)]][['comb.state']][['start']] <- ggplot(df) + geom_line(aes_string(x='lag.in.bp', y='value', col='comb.state', linetype='group')) + geom_text(data=df[seq(1,nrow(df),60),], aes_string(x='lag.in.bp', y='value', label='comb.state'))
		df <- melt(tables.end)
		df$lag.in.bp <- df$lag * binsize
		df$comb.state <- factor(df$comb.state, levels=levels(multi.gr$state))
		ggplts[[as.character(feature)]][['comb.state']][['end']] <- ggplot(df) + geom_line(aes_string(x='lag.in.bp', y='value', col='comb.state')) + geom_text(data=df[seq(1,nrow(df),10),], aes_string(x='lag.in.bp', y='value', label='comb.state'))
		# Tracks
		df <- melt(tracks.start)
		df$lag.in.bp <- df$lag * binsize
		ggplts[[as.character(feature)]][['tracks']][['start']] <- ggplot(df) + geom_line(aes_string(x='lag.in.bp', y='value', col='track'))
		df <- melt(tracks.end)
		df$lag.in.bp <- df$lag * binsize
		ggplts[[as.character(feature)]][['tracks']][['end']] <- ggplot(df) + geom_line(aes_string(x='lag.in.bp', y='value', col='track'))

	}
	

}	
