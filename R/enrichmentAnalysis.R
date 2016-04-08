#' Plot read counts around annotation
#'
#' Plot read counts around annotation as heatmap.
#'
#' @inheritParams enrichmentAtAnnotation
#' @param max.rows An integer specifying the number of randomly subsampled rows that are plotted from the \code{annotation} object. This is necessary to avoid crashing for heatmaps with too many rows.
#' @author Aaron Taudt
#' @importFrom reshape2 melt
#' @export
plotHeatmap <- function(hmm, annotation, bp.around.annotation=10000, max.rows=1000) {

  # Variables
	hmm <- loadMultiHmmsFromFiles(hmm)[[1]]
  binsize <- width(hmm$bins)[1]
  around <- round(bp.around.annotation/binsize)
  
  # Subsampling for plotting of huge data.frames
	annotation <- subsetByOverlaps(annotation, hmm$bins)
  if (length(annotation)>max.rows) {
    annotation <- sample(annotation, size=max.rows, replace=FALSE)
  }
  
	# Get bins that overlap the start of annotation
  ptm <- startTimedMessage("Overlaps with annotation ...")
	index.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], hmm$bins, select="first")
	index.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="last")
	index.plus <- index.plus[!is.na(index.plus)]
	index.minus <- index.minus[!is.na(index.minus)]
	index <- c(index.plus, index.minus)
	stopTimedMessage(ptm)
	
	# Get surrounding indices
	ptm <- startTimedMessage("Getting surrounding indices ...")
	ext.index.plus <- array(NA, dim=c(length(index.plus), 2*around+1), dimnames=list(anno=1:length(index.plus), position=binsize*seq(-around, around, 1)))
	for (i1 in 1:length(index.plus)) {
	  ext.index.plus[i1,] <- seq(from=-around+index.plus[i1], to=index.plus[i1]+around)
	}
	if (length(index.minus)>0) {
  	ext.index.minus <- array(NA, dim=c(length(index.minus), 2*around+1), dimnames=list(anno=1:length(index.minus), position=binsize*seq(-around, around, 1)))
  	for (i1 in 1:length(index.minus)) {
  	  ext.index.minus[i1,] <- rev( seq(from=-around+index.minus[i1], to=index.minus[i1]+around) )
  	}
  	ext.index <- rbind(ext.index.plus, ext.index.minus)
	} else {
	  ext.index <- ext.index.plus
	}
	ext.index[ext.index <= 0] <- NA
	rownames(ext.index) <- 1:nrow(ext.index)
	stopTimedMessage(ptm)
	
	## Go through combinations and then tracks to get the read counts
	ptm <- startTimedMessage("Getting read counts")
	counts <- list()
	for (combination in levels(hmm$bins$combination)) {
	  counts[[combination]] <- list()
	  index.combination <- which(hmm$bins$combination[index]==combination)
	  ext.index.combination <- ext.index[index.combination,]
	  if (is.null(dim(ext.index.combination))) {
	    ext.index.combination <- array(ext.index.combination, dim=c(1,dim(ext.index)[[2]]), dimnames=list(anno=rownames(ext.index)[index.combination], position=dimnames(ext.index)[[2]]))
	  }
	  for (ntrack in colnames(hmm$bins$counts)) {
      counts[[combination]][[ntrack]] <- array(hmm$bins$counts[ext.index.combination,ntrack], dim=dim(ext.index.combination), dimnames=dimnames(ext.index.combination))
	  }
	}
	stopTimedMessage(ptm)
	
	## Theme
	custom_theme <- theme(
	  panel.grid = element_blank(),
	  panel.border = element_rect(fill='NA'),
	  panel.background = element_rect(fill='white'),
	  axis.text.y = element_blank()
	)
	
	## Prepare data.frame
	ptm <- startTimedMessage("Making the plot ...")
	# Exclude rare combinations for plotting
	num.comb <- sapply(counts, function(x) { nrow(x[[1]]) })
	comb2keep <- names(num.comb)[num.comb/sum(num.comb) > 0.005]
	counts <- counts[comb2keep]
	df <- reshape2::melt(counts)
	names(df) <- c('id','position','counts','track','combination')
	df$id <- factor(df$id, levels=rev(unique(df$id)))
	df$combination <- factor(df$combination, levels=unique(df$combination))
	df$track <- factor(df$track, levels=colnames(hmm$bins$counts))
	
	## Plot as heatmap
	ggplt <- ggplot(df) + geom_tile(aes_string(x='position', y='id', color='combination'), size=2)
	ggplt <- ggplt + geom_tile(aes_string(x='position', y='id', fill='counts'), alpha=0.6)
	ggplt <- ggplt + facet_wrap( ~ track, nrow=1) + custom_theme
	ggplt <- ggplt + xlab('distance from annotation in [bp]') + ylab('')
	breaks <- sort(c(0, 10^(0:5), max(df$counts, na.rm = TRUE)))
	ggplt <- ggplt + scale_fill_continuous(trans='log1p', breaks=breaks, labels=breaks, low='white', high='black')
	stopTimedMessage(ptm)
	
	return(ggplt)
}


#' Plot enrichment around annotation
#'
#' Plot fold enrichment of combinatorial states around and inside of annotation.
#'
#' @inheritParams enrichmentAtAnnotation
#' @return A \code{\link{ggplot}} object containing the plot.
#' @author Aaron Taudt
#' @importFrom reshape2 melt
#' @export
plotEnrichment <- function(hmm, annotation, bp.around.annotation=10000, region=c("start","inside","end"), num.intervals=20) {

  ### Get fold enrichment
	enrich <- enrichmentAtAnnotation(hmm, annotation, bp.around.annotation=bp.around.annotation, region=region, what='combinations', num.intervals=num.intervals)
	### Prepare for plotting
	df <- reshape2::melt(enrich[[1]])
	df$L1 <- factor(df$L1, levels=c('start','inside','end'))
	df <- rbind(df[df$L1 == 'start',], df[df$L1 == 'inside',], df[df$L1 == 'end',])
	if (length(region)>=2) {
	  df <- df[!(df$L1 == 'start' & df$lag > 0),]
	  df <- df[!(df$L1 == 'end' & df$lag < 0),]
    df$position <- apply(data.frame(df$interval, df$lag), 1, max, na.rm = TRUE)
	} else if (length(region)==1) {
	  df <- df[df$L1 == region,]
	  df$position <- df$lag
	}
  df$position[df$L1 == 'end'] <- df$position[df$L1 == 'end'] + bp.around.annotation
  df$position[df$L1 == 'inside'] <- df$position[df$L1 == 'inside'] * bp.around.annotation

  ### Plot
	ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='combination'), size=2)
	ggplt <- ggplt + theme_bw() + xlab('distance from annotation in [bp]') + ylab('fold enrichment')
	if (length(region)>=2) {
		breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * bp.around.annotation)
		labels <- c(-bp.around.annotation, -bp.around.annotation/2, '0%', '50%', '100%', bp.around.annotation/2, bp.around.annotation)
		ggplt <- ggplt + scale_x_continuous(breaks=breaks, labels=labels)
	}
  return(ggplt)
	
}


#' Enrichment of (combinatorial) states for genomic annotations
#'
#' The function calculates the enrichment of a genomic feature with peaks or combinatorial states. Input is a \code{\link{multiHMM}} object (containing the peak calls and combinatorial states) and a \code{\link{GRanges}} object containing the annotation of interest (e.g. transcription start sites or genes).
#'
#' @author Aaron Taudt
#' @param hmm A \code{\link{multiHMM}} object.
#' @param annotation A \code{\link{GRanges}} object with the annotation of interest.
#' @param bp.around.annotation An integer specifying the number of basepairs up- and downstream of the annotation for which the enrichment will be calculated.
#' @param region A combination of \code{c('start','inside','end')} specifying the region of the annotation for which the enrichment will be calculated. Select \code{'start'} if you have a point-sized annotation like transcription start sites. Select \code{c('start','inside','end')} if you have long annotations like genes.
#' @param what A combination of \code{c('combinations','binstates','counts')} specifying which statistic to calculate.
#' @param num.intervals Number of intervals for enrichment 'inside' of annotation.
#' @return A \code{list()} containing \code{data.frame()}s for enrichment of combinatorial states and binary states at the start, end and inside of the annotation.
#' @importFrom S4Vectors as.factor subjectHits queryHits
enrichmentAtAnnotation <- function(hmm, annotation, bp.around.annotation=10000, region=c('start','inside','end'), what=c('combinations','binstates','counts'), num.intervals=21) {

	if (class(hmm)!=class.multivariate.hmm) {
		if (class(hmm)!='character' | length(hmm)>1) {
			stop("argument 'hmm' expects a multivariate HMM object or a file that contains such an object")
		}
	}
	hmm <- loadMultiHmmsFromFiles(hmm)[[1]]

	## Variables
	binsize <- width(hmm$bins)[1]
	lag <- round(bp.around.annotation/binsize)
	enrich <- list()
	enrich$combinations <- list()
	enrich$binstates <- list()
	enrich$counts <- list()

	## Get combinatorial and binary states
	combinations <- hmm$bins$combination
	tcombinations <- table(combinations)
	if ('binstates' %in% what) {
		binstates <- dec2bin(hmm$bins$state, colnames=hmm$IDs)
	}

	### Calculationg enrichment inside of annotation ###
	if ('inside' %in% region) {
		message("Enrichment inside of annotations")

		intervals <- seq(from=0, to=1, length.out=num.intervals+1)
		widths.annotation <- width(annotation) - 1
		annotation.1bp <- resize(annotation, 1, fix='start')
		# Initialize arrays
		if ('binstates' %in% what) binstates.inside <- array(dim=c(num.intervals+1, length(hmm$IDs)), dimnames=list(interval=intervals, track=hmm$IDs))
		if ('combinations' %in% what) combinations.inside <- array(dim=c(num.intervals+1, length(levels(hmm$bins$combination))), dimnames=list(interval=intervals, combination=levels(hmm$bins$combination)))
		if ('counts' %in% what) counts.inside <- array(dim=c(num.intervals+1, length(hmm$IDs)), dimnames=list(interval=intervals, track=hmm$IDs))

		for (interval in intervals) {
		  shift <- widths.annotation * interval * c(1,-1,1)[as.integer(strand(annotation))]
			shifted.starts <- start(annotation.1bp) + shift
			annotation.shifted <- GRanges(seqnames = seqnames(annotation.1bp), ranges = IRanges(start = shifted.starts, end = shifted.starts), strand = strand(annotation.1bp))
			# Get bins that overlap the shifted annotation
			index.inside.plus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='+' | strand(annotation.shifted)=='*'], hmm$bins, select="first")
			index.inside.minus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='-'], hmm$bins, select="last")
			index.inside.plus <- index.inside.plus[!is.na(index.inside.plus)]
			index.inside.minus <- index.inside.minus[!is.na(index.inside.minus)]
			index <- c(index.inside.plus, index.inside.minus)
			index <- index[index>0 & index<=length(hmm$bins)]
			if ('binstates' %in% what) binstates.inside[as.character(interval),] <- colMeans(binstates[index,])
			if ('combinations' %in% what) {
				fold <- table(combinations[index]) / tcombinations / length(annotation) * length(hmm$bins) # fold enrichment
# 				fold <- table(combinations[index]) / length(annotation) # percentage enrichment
				fold[is.na(fold)] <- 0
				combinations.inside[as.character(interval),] <- fold
			}
			if ('counts' %in% what) counts.inside[as.character(interval),] <- colMeans(hmm$bins$counts[index,])
		}
		if ('binstates' %in% what) {
			enrich$binstates$inside <- binstates.inside
		}
		if ('combinations' %in% what) {
			enrich$combinations$inside <- combinations.inside
		}
		if ('counts' %in% what) {
			enrich$counts$inside <- counts.inside
		}
	}

	### 10000 bp before annotation ###
	if ('start' %in% region) {
		message("Enrichment ",bp.around.annotation,"bp before annotations")
		# Get bins that overlap the start of annotation
		index.start.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], hmm$bins, select="first")
		index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="last")
		index.start.plus <- index.start.plus[!is.na(index.start.plus)]
		index.start.minus <- index.start.minus[!is.na(index.start.minus)]
		# Occurrences at every bin position relative to feature
		if ('binstates' %in% what) binstates.start <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		if ('combinations' %in% what) combinations.start <- array(dim=c(length(-lag:lag), length(levels(hmm$bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(hmm$bins$combination)))
		if ('counts' %in% what) counts.start <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		for (ilag in -lag:lag) {
			index <- c(index.start.plus+ilag, index.start.minus-ilag)
			index <- index[index>0 & index<=length(hmm$bins)]
			if ('binstates' %in% what) binstates.start[as.character(ilag),] <- colMeans(binstates[index,])
			if ('combinations' %in% what) {
				fold <- table(combinations[index]) / tcombinations / length(annotation) * length(hmm$bins) # fold enrichment
# 				fold <- table(combinations[index]) / length(annotation) # percentage enrichment
				fold[is.na(fold)] <- 0
				combinations.start[as.character(ilag),] <- fold
			}
			if ('counts' %in% what) counts.start[as.character(ilag),] <- colMeans(hmm$bins$counts[index,])
		}
		if ('binstates' %in% what) {
			rownames(binstates.start) <- as.numeric(rownames(binstates.start)) * binsize
			enrich$binstates$start <- binstates.start
		}
		if ('combinations' %in% what) {
			rownames(combinations.start) <- as.numeric(rownames(combinations.start)) * binsize
			enrich$combinations$start <- combinations.start
		}
		if ('counts' %in% what) {
			rownames(counts.start) <- as.numeric(rownames(counts.start)) * binsize
			enrich$counts$start <- counts.start
		}
	}

	### 10000 bp after annotation ###
	if ('end' %in% region) {
		message("Enrichment ",bp.around.annotation,"bp after annotations")
		# Get bins that overlap the end of annotation
		index.end.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], hmm$bins, select="last")
		index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], hmm$bins, select="first")
		index.end.plus <- index.end.plus[!is.na(index.end.plus)]
		index.end.minus <- index.end.minus[!is.na(index.end.minus)]
		# Occurrences at every bin position relative to feature
		if ('binstates' %in% what) binstates.end <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		if ('combinations' %in% what) combinations.end <- array(dim=c(length(-lag:lag), length(levels(hmm$bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(hmm$bins$combination)))
		if ('counts' %in% what) counts.end <- array(dim=c(length(-lag:lag), length(hmm$IDs)), dimnames=list(lag=-lag:lag, track=hmm$IDs))
		for (ilag in -lag:lag) {
	# 		message("lag = ",ilag)
			index <- c(index.end.plus+ilag, index.end.minus-ilag)
			index <- index[index>0 & index<=length(hmm$bins)]
			if ('binstates' %in% what) binstates.end[as.character(ilag),] <- colMeans(binstates[index,])
			if ('combinations' %in% what) {
				fold <- table(combinations[index]) / tcombinations / length(annotation) * length(hmm$bins) # fold enrichment
# 				fold <- table(combinations[index]) / length(annotation) # percentage enrichment
				fold[is.na(fold)] <- 0
				combinations.end[as.character(ilag),] <- fold
			}
			if ('counts' %in% what) counts.end[as.character(ilag),] <- colMeans(hmm$bins$counts[index,])
		}
		if ('binstates' %in% what) {
			rownames(binstates.end) <- as.numeric(rownames(binstates.end)) * binsize
			enrich$binstates$end <- binstates.end
		}
		if ('combinations' %in% what) {
			rownames(combinations.end) <- as.numeric(rownames(combinations.end)) * binsize
			enrich$combinations$end <- combinations.end
		}
		if ('counts' %in% what) {
			rownames(counts.end) <- as.numeric(rownames(counts.end)) * binsize
			enrich$counts$end <- counts.end
		}
	}

	return(enrich)

}

