plot.cross.correlation <- function(hmm, annotation, bp.around.annotation=10000) {

	## Debugging
	library(biomaRt)
	library(GenomicRanges)
	library(ggplot2)
	library(reshape2)
	library(chromstar)
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

	## Libraries
	library(ggplot2)
	library(reshape2)

	## Variables
	binsize <- width(hmm$bins)[1]
	lag <- round(bp.around.annotation/binsize)

	### Calculationg cross-correlation inside of normed annotation ###
	cat("calculating cross-correlation inside of annotation\n")
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
	cat("calculating cross-correlation before annotation\n")
	# Get bins that overlap the start of annotation
	index.start.plus <- findOverlaps(annotation[strand(annotation)=='+'], hmm$bins, select="first")
	index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="last")
	index.start.plus <- index.start.plus[!is.na(index.start.plus)]
	index.start.minus <- index.start.minus[!is.na(index.start.minus)]
	# Occurrences at every bin position relative to feature
	binstates.before <- array(dim=c(length(-lag:0), ncol(hmm$bins$reads)), dimnames=list(lag=-lag:0, track=colnames(hmm$bins$reads)))
	reads.before <- array(dim=c(length(-lag:0), ncol(hmm$bins$reads)), dimnames=list(lag=-lag:0, track=colnames(hmm$bins$reads)))
	for (ilag in -lag:0) {
		cat("lag =",ilag,"\r")
		index <- c(index.start.plus+ilag, index.start.minus-ilag)
		index <- index[index>0]
		binstates.before[as.character(ilag),] <- colMeans(binstates[index,])
		reads.before[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		cat("             \r")
	}
	cat("calculating cross-correlation after annotation\n")
	# Get bins that overlap the end of annotation
	index.end.plus <- findOverlaps(annotation[strand(annotation)=='+'], hmm$bins, select="last")
	index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="first")
	index.end.plus <- index.end.plus[!is.na(index.end.plus)]
	index.end.minus <- index.end.minus[!is.na(index.end.minus)]
	# Occurrences at every bin position relative to feature
	binstates.after <- array(dim=c(length(0:lag), ncol(hmm$bins$reads)), dimnames=list(lag=0:lag, track=colnames(hmm$bins$reads)))
	reads.after <- array(dim=c(length(0:lag), ncol(hmm$bins$reads)), dimnames=list(lag=0:lag, track=colnames(hmm$bins$reads)))
	for (ilag in 0:lag) {
		cat("lag =",ilag,"\r")
		index <- c(index.end.plus+ilag, index.end.minus-ilag)
		index <- index[index>0]
		binstates.after[as.character(ilag),] <- colMeans(binstates[index,])
		reads.after[as.character(ilag),] <- colMeans(hmm$bins$reads[index,])
		cat("             \r")
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
	ggplt.reads <- ggplot(df) + geom_line(aes(x=positions.plot, y=value, col=track))
	df <- melt(binstates.around, varnames=c('pos','track'), id.vars=c('positions','positions.plot'), variable.name='track')
	ggplt.binstates <- ggplot(df) + geom_line(aes(x=positions.plot, y=value, col=track))


	
}

cross.correlation <- function(multi.hmm, annotation.file, grouping=NULL, lag.in.bp=10000) {

	## Import annotation file
	cat("importing",annotation.file," ...")
	anno.data <- rtracklayer::import(annotation.file)
	cat(" done\n")

	## Convert multi.hmm to GRanges
	if (class(multi.hmm) == class.multivariate.hmm) {
		cat("converting multi.hmm to GRanges...")
		multi.gr <- hmm2GRanges(multi.hmm, reduce=F)
		cat(" done\n")
	} else if (class(multi.hmm) == 'GRanges') {
		multi.gr <- multi.hmm
	} else {
		stop("argument 'multi.hmm' expects multivariate HMM object")
	}

	## Calculate variables
	features <- levels(anno.data$type)
	binsize <- width(multi.gr[1])
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
		cat(paste0("Calculating cross correlation for feature '",feature,"'\n"))
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
			cat("lag =",ilag,"\r")
			# Combinatorial states
			for (group in groups) {
				tables.start[as.character(ilag), ,as.character(group)] <- table(states.per.group[[as.character(group)]][index.start.plus+ilag]) + table(states.per.group[[as.character(group)]][index.start.minus-ilag])
				tables.end[as.character(ilag), ,as.character(group)] <- table(states.per.group[[as.character(group)]][index.end.plus+ilag]) + table(states.per.group[[as.character(group)]][index.end.minus-ilag])
			}
			# Tracks
			tracks.start[as.character(ilag),] <- colSums(dec2bin(multi.gr[index.start.plus+ilag]$state, ndigits=ncol(multi.gr$reads))) + colSums(dec2bin(multi.gr[index.start.minus-ilag]$state, ndigits=ncol(multi.gr$reads)))
			tracks.end[as.character(ilag),] <- colSums(dec2bin(multi.gr[index.end.plus+ilag]$state, ndigits=ncol(multi.gr$reads))) + colSums(dec2bin(multi.gr[index.end.minus-ilag]$state, ndigits=ncol(multi.gr$reads)))
			cat("             \r")
		}

		## Plot cross-correlations
		library(ggplot2)
		library(reshape2)
		# Combinatorial states
		df <- melt(tables.start)
		df$lag.in.bp <- df$lag * binsize
		df$comb.state <- factor(df$comb.state, levels=levels(multi.gr$state))
		df$group <- factor(df$group, levels=levels.group)
		ggplts[[as.character(feature)]][['comb.state']][['start']] <- ggplot(df) + geom_line(aes(x=lag.in.bp, y=value, col=comb.state, linetype=group)) + geom_text(data=df[seq(1,nrow(df),60),], aes(x=lag.in.bp, y=value, label=comb.state))
		df <- melt(tables.end)
		df$lag.in.bp <- df$lag * binsize
		df$comb.state <- factor(df$comb.state, levels=levels(multi.gr$state))
		ggplts[[as.character(feature)]][['comb.state']][['end']] <- ggplot(df) + geom_line(aes(x=lag.in.bp, y=value, col=comb.state)) + geom_text(data=df[seq(1,nrow(df),10),], aes(x=lag.in.bp, y=value, label=comb.state))
		# Tracks
		df <- melt(tracks.start)
		df$lag.in.bp <- df$lag * binsize
		ggplts[[as.character(feature)]][['tracks']][['start']] <- ggplot(df) + geom_line(aes(x=lag.in.bp, y=value, col=track))
		df <- melt(tracks.end)
		df$lag.in.bp <- df$lag * binsize
		ggplts[[as.character(feature)]][['tracks']][['end']] <- ggplot(df) + geom_line(aes(x=lag.in.bp, y=value, col=track))

	}
	

}	
