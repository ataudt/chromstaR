#' Combine combinatorial states from several Multivariates
#' 
#' Combine combinatorial states from several \code{\link{multiHMM}} objects. Combinatorial states can be combined for objects containing multiple marks (\code{mode='mark'}) or multiple conditions (\code{mode='condition'}).
#' 
#' @param multi.hmm.list A named \code{list()} with \code{\link{multiHMM}} objects. The names of the list are used to name the conditions or marks. Alternatively a named character vector with filenames that contain \code{\link{multiHMM}} objects.
#' @param mode Mode of combination. See \code{\link{Chromstar}} for a description of the \code{mode} parameter.
#' @param conditions A character vector with conditions in case of \code{mode='full'} or \code{mode='condition'}.
#' @return A \code{link{combinedHMM}} objects with combinatorial states for each condition.
#' @author Aaron Taudt
#' @export
combineMultivariates <- function(multi.hmm.list, mode, conditions) {
	
  ## Helper function
  getCondition <- function(combinations, condition) {
		combinations <- sub('\\[','',combinations)
		combinations <- sub('\\]','',combinations)
	  combinations.split <- strsplit(as.character(combinations), '\\+')
    condition.string <- paste0('\\<', condition, '\\>')
    combinations.split.condition <- lapply(combinations.split, function(x) { grep(condition.string, x, value=TRUE) })
    combinations.condition <- sapply(combinations.split.condition, paste, collapse='+')
    combinations.condition <- gsub(paste0('-',condition.string), '', combinations.condition)
    combinations.condition <- paste0('[', combinations.condition, ']')
    return(combinations.condition)
  }
	  
	if (mode == 'mark') {
		if (is.null(names(multi.hmm.list))) {
			stop("'multi.hmm.list' must be named.")
		}
		## Load first HMM for coordinates
		ptm <- startTimedMessage("Processing condition ",names(multi.hmm.list)[1]," ...")
		hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[1]])[[1]] )
		bins <- hmm$bins
		mcols(bins) <- NULL
		## Add combinatorial states
		combs <- list()
		combs[[names(multi.hmm.list)[1]]] <- hmm$bins$combination
		stopTimedMessage(ptm)
		
		for (i1 in 2:length(multi.hmm.list)) {
			ptm <- startTimedMessage("Processing condition ",names(multi.hmm.list)[i1]," ...")
			hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[i1]])[[1]] )
			combs[[names(multi.hmm.list)[i1]]] <- hmm$bins$combination
			stopTimedMessage(ptm)
		}
		combs.df <- as(combs,'DataFrame')
		
	} else if (mode == 'condition') {
		if (is.null(names(multi.hmm.list))) {
			stop("'multi.hmm.list' must be named.")
		}
		### Get vectors with presence/absence of each mark and condition
		states <- list()
		for (mark in names(multi.hmm.list)) {
			ptm <- startTimedMessage("Processing mark ",mark," ...")
			i1 <- which(mark==names(multi.hmm.list))
			## Load HMM
			hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[mark]])[[1]] )
			## Extract conditions
			comblevels <- levels(hmm$bins$combination)
			comblevels <- sub('\\[','', comblevels)
			comblevels <- sub('\\]','', comblevels)
			conds.split <- strsplit(as.character(comblevels),'\\+')
			conds <- conds.split[[which.max(sapply(conds.split, length))]]
			if (i1==1) {
				bins <- hmm$bins
				mcols(bins) <- NULL
			}
			if (any(!conds %in% conditions)) {
				stop(paste0("Not all levels(multi.hmm.list[[",i1,"]]$bins$combination) in 'conditions'."))
			}
			## Make comblevel -> mark mappings and map differential states to mark
			states[[mark]] <- list()
			for (cond in conditions) {
				mapping <- grepl(paste0('\\<',cond,'\\>'), conds.split)
				names(mapping) <- comblevels
				states[[mark]][[cond]] <- c('',mark)[mapping[hmm$bins$combination]+1] # no need to coerce to character here because the order is the same
			}
			stopTimedMessage(ptm)
		}
		### Paste the marks over each condition
		message("Pasting into combinatorial states")
		combs <- list()
		for (cond in conditions) {
			ptm <- startTimedMessage("  condition ",cond," ...")
			l <- lapply(states, '[[', cond)
			l$sep <- '+'
			comb <- do.call(paste, l)
			comb <- gsub('\\+{2,}','+', comb)
			comb <- sub('^\\+','', comb)
			comb <- sub('\\+$','', comb)
			comb <- paste0('[', comb, ']')
			combs[[cond]] <- comb
			stopTimedMessage(ptm)
		}
		combs.df <- as.data.frame(combs)
		combs.df <- as(combs.df, 'DataFrame')

	} else if (mode == 'full') {
	  if (length(multi.hmm.list) > 1) {
	    stop("'multi.hmm.list' must contain only one 'multiHMM' object when mode='full'.")
	  }
	  
		hmm <- suppressMessages( loadMultiHmmsFromFiles(multi.hmm.list[[1]])[[1]] )
		bins <- hmm$bins
		mcols(bins) <- NULL
	  combs <- list()
	  for (condition in conditions) {
  		ptm <- startTimedMessage("Processing condition ",condition," ...")
      mapping.condition <- getCondition(hmm$mapping, condition)
	    names(mapping.condition) <- names(hmm$mapping)
	    
	    combs[[as.character(condition)]] <- mapping.condition[hmm$bins$state]
	    # # Assign levels
	    # comblevels.condition <- unique(getCondition(levels(hmm$bins$combination), condition))
	    # levels(combs[[as.character(condition)]]) <- comblevels.condition
	    stopTimedMessage(ptm)
	  }
	  combs.df <- as.data.frame(combs)
	  combs.df <- as(combs.df, 'DataFrame')
		
	} else {
		stop("Please specify either 'conditions' or 'marks'.")
	}
  # Reassign levels such that all conditions have the same levels
	comblevels <- sort(unique(unlist(lapply(combs.df, levels))))
	combs.df <- endoapply(combs.df, function(x) { x <- factor(x, levels=comblevels) })
	
	### Redo the segmentation for all conditions combined
	ptm <- startTimedMessage("Redoing segmentation for all conditions combined ...")
	values(bins) <- combs.df
	l <- as.list(combs.df)
	mcols(bins)$state <- as.numeric(do.call(paste0, lapply(l, as.integer)))
	bins.df <- as.data.frame(bins)
	segments.df <- suppressMessages( collapseBins(bins.df, column2collapseBy='state', columns2drop=c('width','state')) )
	combined.segments <- as(segments.df, 'GRanges')
	seqlengths(combined.segments) <- seqlengths(bins)
	mcols(bins)$state <- NULL
	stopTimedMessage(ptm)
	
	### Redo the segmentation for each condition separately
	ptm <- startTimedMessage("Redoing segmentation for each condition separately ...")
	segments <- list()
	for (cond in names(combs)) {
		bins.cond <- bins
		mcols(bins.cond) <- mcols(bins)[cond]
		df <- as.data.frame(bins.cond)
		segments.cond <- suppressMessages( collapseBins(df, column2collapseBy=cond, columns2drop=c('width')) )
		segments.cond <- as(segments.cond, 'GRanges')
		names(mcols(segments.cond)) <- 'combination'
		seqlengths(segments.cond) <- seqlengths(bins)
		segments[[cond]] <- segments.cond
	}
	stopTimedMessage(ptm)
	
	### Make return object
	hmm <- list()
	class(hmm) <- class.combined.multivariate.hmm
	hmm$bins <- bins
	hmm$combined.segments <- combined.segments
	hmm$segments <- segments
	return(hmm)
	
}
