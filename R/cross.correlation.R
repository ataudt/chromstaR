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
