#' Combine combinatorial states from several Multivariates
#' 
#' Combine combinatorial states from several \code{\link{chromstaR_multivariateHMM}} objects. Combinatorial states can be combined across \code{conditions} (from a combinatorial analysis) or across \code{marks} (from a differential analysis).
#' 
#' @param conditions A named \code{list()} with \code{\link{chromstaR_multivariateHMM}} objects. The names of the list are used to name the conditions. Alternatively a named character vector with filenames that contain \code{\link{chromstaR_multivariateHMM}} objects.
#' @param marks A named \code{list()} with \code{\link{chromstaR_multivariateHMM}} objects. The names of the list are used to name the marks. Alternatively a named character vector with filenames that contain \code{\link{chromstaR_multivariateHMM}} objects.
#' @return A DataFrame with combinatorial states for each condition.
#' @author Aaron Taudt
#' @export
combineMultivariates <- function(conditions=NULL, marks=NULL) {
	
	if (!is.null(conditions) & !is.null(marks)) {
		stop("Both 'conditions' and 'marks' are specified. Please specify only one of those.")
	}
	if (length(conditions)<=1 & !is.null(conditions)) {
		stop("'conditions' must contain at least two entries.")
	}
	if (length(marks)<=1 & !is.null(marks)) {
		stop("'marks' must contain at least two entries.")
	}
	
	if (!is.null(conditions)) {
		if (is.null(names(conditions))) {
			stop("'conditions' must be named.")
		}
		## Load first HMM for coordinates
		message("Processing condition ",names(conditions)[1]," ...", appendLF=FALSE); ptm <- proc.time()
		hmm <- suppressMessages( loadMultiHmmsFromFiles(conditions[[1]])[[1]] )
		bins <- hmm$bins
		mcols(bins) <- NULL
		## Add combinatorial states
		combs <- list()
		combs[[names(conditions)[1]]] <- hmm$bins$combination
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		
		for (i1 in 2:length(conditions)) {
			message("Processing condition ",names(conditions)[i1]," ...", appendLF=FALSE); ptm <- proc.time()
			hmm <- suppressMessages( loadMultiHmmsFromFiles(conditions[[i1]])[[1]] )
			combs[[names(conditions)[i1]]] <- hmm$bins$combination
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		combs.df <- as(combs,'DataFrame')
		
	} else if (!is.null(marks)) {
		if (is.null(names(marks))) {
			stop("'marks' must be named.")
		}
		### Get vectors with presence/absence of each mark and condition
		states <- list()
		for (mark in names(marks)) {
			message("Processing mark ",mark," ...", appendLF=FALSE); ptm <- proc.time()
			i1 <- which(mark==names(marks))
			## Load HMM
			hmm <- suppressMessages( loadMultiHmmsFromFiles(marks[[mark]])[[1]] )
			## Extract conditions
			comblevels <- levels(hmm$bins$combination)
			conds.split <- strsplit(comblevels,'-')
			conds <- conds.split[[which.max(sapply(conds.split, length))]]
			if (i1==1) {
				bins <- hmm$bins
				mcols(bins) <- NULL
				comblevels.1 <- comblevels
				conds.1 <- conds
			}
			if (!all.equal(conds.1, conds)) {
				stop(paste0("levels(marks[[x]]$bins$combination) differ between x=1 and x=",i1,". They should be the same. Please check your input."))
			}
			## Make comblevel -> mark mappings and map differential states to mark
			states[[mark]] <- list()
			for (cond in conds) {
				mapping <- grepl(paste0('\\<',cond,'\\>'), conds.split)
				names(mapping) <- comblevels
				states[[mark]][[cond]] <- c('',mark)[mapping[hmm$bins$combination]+1] # no need to coerce to character here because the order is the same
			}
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		### Paste the marks over each condition
		message("Pasting into combinatorial states")
		combs <- list()
		for (cond in conds) {
			message("  condition ",cond," ...", appendLF=FALSE); ptm <- proc.time()
			l <- lapply(states, '[[', cond)
			l$sep <- '+'
			comb <- do.call(paste, l)
			comb <- gsub('\\+{2,}','+', comb)
			comb <- sub('^\\+','', comb)
			comb <- sub('\\+$','', comb)
			combs[[cond]] <- comb
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		combs.df <- as.data.frame(combs)
		combs.df <- as(combs.df, 'DataFrame')

	} else {
		stop("Please specify either 'conditions' or 'marks'.")
	}
	
	### Redo the segmentation for all conditions combined
	message("Redoing segmentation for all conditions combined ...", appendLF=FALSE); ptm <- proc.time()
	values(bins) <- combs.df
	l <- as.list(combs.df)
	mcols(bins)$state <- as.numeric(do.call(paste0, lapply(l, as.integer)))
	bins.df <- as.data.frame(bins)
	segments.df <- suppressMessages( collapseBins(bins.df, column2collapseBy='state', columns2drop=c('width','state')) )
	combined.segments <- as(segments.df, 'GRanges')
	seqlengths(combined.segments) <- seqlengths(bins)
	mcols(bins)$state <- NULL
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	
	### Redo the segmentation for each condition separately
	message("Redoing segmentation for each condition separately ...", appendLF=FALSE); ptm <- proc.time()
	segments <- list()
	for (cond in names(combs)) {
		bins.cond <- bins
		mcols(bins.cond) <- mcols(bins)[cond]
		df <- as.data.frame(bins.cond)
		segments.cond <- suppressMessages( collapseBins(df, column2collapseBy=cond, columns2drop=c('width')) )
		segments.cond <- as(segments.cond, 'GRanges')
		seqlengths(segments.cond) <- seqlengths(bins)
		segments[[cond]] <- segments.cond
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	
	### Make return object
	hmm <- list()
	class(hmm) <- class.combined.multivariate.hmm
	hmm$bins <- bins
	hmm$combined.segments <- combined.segments
	hmm$segments <- segments
	
}