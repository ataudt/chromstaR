#' Fold enrichment of combinatorial states
#'
#' Compute the fold enrichment of combinatorial states in a given feature (e.g. TSS, exons, ...)
#'
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} or a file that contains such an object.
#' @param feature A \code{\link{GRanges}} with coordinates of the feature to compute the fold enrichment to.
#' @param featurelist A list with \code{\link{Granges}} objects containing coordinates of multiple features. The names of the list entries will be used to name the return values.
#' @param combinations A vector with combinations for which the fold enrichment will be calculated. If \code{NULL} all combinations will be considered.
#' @param percentages Set to \code{TRUE} if you want to have percentages (0 to 1) instead of fold enrichments returned. Note that in this case different features are not directly comparable.
#' @return A named array with fold enrichments. If \code{percentages=TRUE} a list with arrays of percentage (0 to 1) enrichments.
#' @author Aaron Taudt
#' @export
foldEnrichment <- function(multi.hmm, featurelist, combinations=NULL, percentages=FALSE) {
	
	multi.hmm <- loadMultiHmmsFromFiles(multi.hmm)[[1]]
	## Variables
	bins <- multi.hmm$bins
	if (is.null(combinations)) {
		comb.levels <- levels(bins$combination)
	} else {
		comb.levels <- combinations
	}
	genome <- sum(as.numeric(width(bins)))
	feature.lengths <- lapply(featurelist, function(x) { sum(as.numeric(width(x))) })
	
	## Fold enrichment
	fold <- array(NA, dim=c(length(featurelist), length(comb.levels)), dimnames=list(feature=names(featurelist), combination=comb.levels))
	perc.combstate.in.feature <- fold
	perc.feature.in.combstate <- fold
	for (icomb in 1:length(comb.levels)) {
		mask <- bins$combination == comb.levels[icomb]
		bins.mask <- bins[mask]
		combstate.length <- sum(as.numeric(width(bins.mask)))
		for (ifeat in 1:length(featurelist)) {
			feature <- featurelist[[ifeat]]
			ind <- findOverlaps(bins.mask, feature)

			binsinfeature <- bins.mask[unique(queryHits(ind))]
			sum.binsinfeature <- sum(as.numeric(width(binsinfeature)))
			perc.combstate.in.feature[ifeat,icomb] <- sum.binsinfeature / combstate.length

			featuresinbins <- feature[unique(subjectHits(ind))]
			sum.featuresinbins <- sum(as.numeric(width(featuresinbins)))
			perc.feature.in.combstate[ifeat,icomb] <- sum.featuresinbins / feature.lengths[[ifeat]]

			fold[ifeat,icomb] <- sum.binsinfeature / combstate.length / feature.lengths[[ifeat]] * genome
		}
	}
	
	if (percentages) {
		return(list(combstate.in.feature=perc.combstate.in.feature, feature.in.combstate=perc.feature.in.combstate))
	}
	return(fold)
}


#' Overlap with expression data
#'
#' Get the expression values that overlap with each combinatorial state.
#'
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} or a file that contains such an object.
#' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range.
#' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
#' @param return.marks Set to \code{TRUE} if expression values for marks instead of combinations should be returned.
#' @return A named list with expression values.
#' @author Aaron Taudt
#' @export
expressionOverlap <- function(multi.hmm, expression, combinations=NULL, return.marks=FALSE) {
	
	multi.hmm <- loadMultiHmmsFromFiles(multi.hmm)[[1]]
	## Variables
	bins <- multi.hmm$bins
	if (is.null(combinations)) {
		comb.levels <- levels(bins$combination)
	} else {
		comb.levels <- combinations
	}
	marks <- mixedsort(unique(unlist(strsplit(comb.levels,'-'))))
	
	exprlist <- list()
	if (return.marks) {
		for (mark in marks) {
			mask <- grepl(paste0('\\<',mark,'\\>'),bins$combination)
			expr.mark <- subsetByOverlaps(expression, bins[mask])
			exprlist[[mark]] <- expr.mark$expression
		}
	} else {
		for (comb.level in comb.levels) {
			mask <- bins$combination == comb.level
			expr.combstate <- subsetByOverlaps(expression, bins[mask])
			exprlist[[comb.level]] <- expr.combstate$expression
		}
	}
	return(exprlist)

}


#' Mean expression at percentage overlap
#'
#' Get the average expression for each percentage of overlap of combinatorial state with feature.
#'
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} or a file that contains such an object.
#' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range of the feature.
#' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
#' @return A list with vectors of mean expression values per percentile for each combinatorial state. 
#' @author Aaron Taudt
#' @export
expressionAtPercentageOverlap <- function(multi.hmm, expression, combinations=NULL) {

	multi.hmm <- loadMultiHmmsFromFiles(multi.hmm)[[1]]
	## Variables
	bins <- multi.hmm$bins
	if (is.null(combinations)) {
		comb.levels <- levels(bins$combination)
	} else {
		comb.levels <- combinations
	}
	nintervals <- 100
	
	expression.means <- array(NA, dim=c(nintervals+1, length(comb.levels), 2), dimnames=list(percentage=0:nintervals, combination=comb.levels, value=c('expression','weight')))
	for (icomb in 1:length(comb.levels)) {
		mask <- bins$combination == comb.levels[icomb]
		ind <- findOverlaps(expression, bins[mask])
		rle <- rle(queryHits(ind))
		expression$num.bins <- 0
		expression$num.bins[rle$values] <- rle$lengths
		expression$genewidth <- width(expression)
		expression$percentage <- round(expression$num.bins*1000 / expression$genewidth * nintervals) # Normalize to genewidth
		expression$percentage[expression$percentage>=nintervals] <- nintervals
		splt <- split(expression$expression, expression$percentage)
		tab <- sapply(splt, mean, na.rm=TRUE)
		expression.means[names(tab),icomb,'expression'] <- tab #select by icomb instead of name because of potential '' states
		expression.means[names(tab),icomb,'weight'] <- sapply(splt, length)
	}
	return(expression.means)

}


#' Frequency of combinatorial states
#'
#' Get the genomewide frequency of each combinatorial state.
#'
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} or a file that contains such an object.
#' @param combinations A vector with combinations for which the frequency will be calculated. If \code{NULL} all combinations will be considered.
#' @return A table with frequencies of each combinatorial state.
#' @author Aaron Taudt
#' @export
combinatorialFrequency <- function(multi.hmm, combinations=NULL) {

	multi.hmm <- loadMultiHmmsFromFiles(multi.hmm)[[1]]
	## Variables
	bins <- multi.hmm$bins
	if (is.null(combinations)) {
		comb.levels <- levels(bins$combination)
	} else {
		comb.levels <- combinations
	}

	t <- table(bins$combination) / length(bins)
	t <- t[names(t) %in% comb.levels]

	return(t)
}


#' Transition frequencies of combinatorial states
#'
#' Get a table of transition frequencies between combinatorial states of different \code{\link{chromstaR_multivariateHMM}}s.
#'
#' @param multi.hmms A list with \code{\link{chromstaR_multivariateHMM}} objects or a vector with filenames that contain such objects.
#' @param zero.states The string(s) which identifies the zero.states.
#' @param combstates Alternative input instead of \code{multi.hmms}: A list of combinatorial state vectors instead of HMMs. If this is specified, \code{multi.hmms} will be ignored.
#' @return A data.frame with transition frequencies.
#' @author Aaron Taudt
#' @export
transitionFrequencies <- function(multi.hmms=NULL, zero.states="", combstates=NULL) {

	if (is.null(combstates)) {
		## Get combinatorial states in loop to save memory
		combstates <- list()
		for (imodel in 1:length(multi.hmms)) {
			multi.hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmms[[imodel]])[[1]] )
			combstates[[imodel]] <- multi.hmm$bins$combination
		}
		names(combstates) <- names(multi.hmms)
	}
	num.models <- length(combstates)

	### Get transitions for whole genome ###
	combstates$sep <- '>'
	gentrans <- do.call(paste, combstates)
	freqtrans <- as.data.frame(table(gentrans) / length(gentrans))
	names(freqtrans) <- c('transition','frequency')
	freqtrans <- freqtrans[order(freqtrans$frequency, decreasing=TRUE),]
	# Cumulative frequencies
	freqtrans$cumulative.frequency <- cumsum(freqtrans$frequency)

	### Assigning groups for frequency table ###
	freqtrans$group <- 'other'
	freqtrans$transition <- gsub('-','_',freqtrans$transition) # words separated with underscore will be treated as one when searching with \\<XXX\\>
	# Stage-specific and constant states
	combstates$sep <- NULL
	levels.combstates <- unique(unlist(lapply(combstates, levels)))
	levels.combstates <- setdiff(levels.combstates, zero.states)
	levels.combstates <- gsub('-','_',levels.combstates)
	for (combination in levels.combstates) {
		mask <- intersect(grep(paste0('\\<',combination,'\\>'), freqtrans$transition), grep(paste(paste0('\\<',setdiff(levels.combstates,combination),'\\>'), collapse='|'), freqtrans$transition, invert=T))
		freqtrans$group[mask] <- paste0('stage-specific ',combination)
		mask <- sapply(gregexpr(paste0('\\<',combination,'\\>'),freqtrans$transition), length) == num.models
		freqtrans$group[mask] <- paste0('constant ', combination)
	}
	# Zero transitions
	freqtrans.split <- strsplit(sub('>$','>>',as.character(freqtrans$transition)),'>')
	freqtrans.split <- do.call(rbind, freqtrans.split)
	df <- as.data.frame(apply(freqtrans.split, 2, function(x) { x %in% zero.states } ))
	if (ncol(df)==1) {
		iszero <- Reduce('&', df[,1])
	} else {
		iszero <- Reduce('&', as.list(df))
	}
	freqtrans$group[iszero] <- 'zero transition'

	### Assigning groups over whole genome ###
	mapping <- freqtrans$group
	names(mapping) <- freqtrans$transition
	gengroups <- mapping[gentrans]
	gentrans <- data.frame(transition=gentrans, group=gengroups)


	## Return value ##
	return(list(frequencies=freqtrans, transitions=gentrans))

}
