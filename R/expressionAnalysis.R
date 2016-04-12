# #' Overlap with expression data
# #'
# #' Get the expression values that overlap with each combinatorial state.
# #'
# #' @param multi.hmm A \code{\link{multiHMM}} or a file that contains such an object.
# #' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range.
# #' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
# #' @param return.marks Set to \code{TRUE} if expression values for marks instead of combinations should be returned.
# #' @return A named list with expression values.
# #' @author Aaron Taudt
# #' @importFrom gtools mixedsort
# #' @export
# expressionOverlap <- function(multi.hmm, expression, combinations=NULL, return.marks=FALSE) {
# 	
# 	multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=class.multivariate.hmm)[[1]]
# 	## Variables
# 	bins <- multi.hmm$bins
# 	if (is.null(combinations)) {
# 		comb.levels <- levels(bins$combination)
# 	} else {
# 		comb.levels <- combinations
# 	}
# 	marks <- gtools::mixedsort(unique(unlist(strsplit(comb.levels,'+'))))
# 	
# 	exprlist <- list()
# 	if (return.marks) {
# 		for (mark in marks) {
# 			mask <- grepl(paste0('\\<',mark,'\\>'),bins$combination)
# 			expr.mark <- subsetByOverlaps(expression, bins[mask])
# 			exprlist[[mark]] <- expr.mark$expression
# 		}
# 	} else {
# 		for (comb.level in comb.levels) {
# 			mask <- bins$combination == comb.level
# 			expr.combstate <- subsetByOverlaps(expression, bins[mask])
# 			exprlist[[comb.level]] <- expr.combstate$expression
# 		}
# 	}
# 	return(exprlist)
# 
# }
# 
# 
# #' Mean expression at percentage overlap
# #'
# #' Get the average expression for each percentage of overlap of combinatorial state with feature.
# #'
# #' @param multi.hmm A \code{\link{multiHMM}} or a file that contains such an object.
# #' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range of the feature.
# #' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
# #' @return A list with vectors of mean expression values per percentile for each combinatorial state. 
# #' @author Aaron Taudt
# #' @importFrom S4Vectors queryHits
# #' @export
# expressionAtPercentageOverlap <- function(multi.hmm, expression, combinations=NULL) {
# 
# 	multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=class.multivariate.hmm)[[1]]
# 	## Variables
# 	bins <- multi.hmm$bins
# 	if (is.null(combinations)) {
# 		comb.levels <- levels(bins$combination)
# 	} else {
# 		comb.levels <- combinations
# 	}
# 	nintervals <- 100
# 	
# 	expression.means <- array(NA, dim=c(nintervals+1, length(comb.levels), 2), dimnames=list(percentage=0:nintervals, combination=comb.levels, value=c('expression','weight')))
# 	for (icomb in 1:length(comb.levels)) {
# 		mask <- bins$combination == comb.levels[icomb]
# 		ind <- findOverlaps(expression, bins[mask])
# 		rle <- rle(S4Vectors::queryHits(ind))
# 		expression$num.bins <- 0
# 		expression$num.bins[rle$values] <- rle$lengths
# 		expression$genewidth <- width(expression)
# 		expression$percentage <- round(expression$num.bins*1000 / expression$genewidth * nintervals) # Normalize to genewidth
# 		expression$percentage[expression$percentage>=nintervals] <- nintervals
# 		splt <- split(expression$expression, expression$percentage)
# 		tab <- sapply(splt, mean, na.rm=TRUE)
# 		expression.means[names(tab),icomb,'expression'] <- tab #select by icomb instead of name because of potential '' states
# 		expression.means[names(tab),icomb,'weight'] <- sapply(splt, length)
# 	}
# 	return(expression.means)
# 
# }
