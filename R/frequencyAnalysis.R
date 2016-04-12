#' Frequencies of combinatorial states
#'
#' Get the genomewide frequency of each combinatorial state.
#'
#' @param multi.hmm A \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object or a file that contains such an object.
#' @param combinations A vector with combinations for which the frequency will be calculated. If \code{NULL} all combinations will be considered.
#' @return A table with frequencies of each combinatorial state.
#' @author Aaron Taudt
#' @export
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-mark_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'genomicFrequencies(model)
#'
genomicFrequencies <- function(multi.hmm, combinations=NULL) {

	multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
	bins <- multi.hmm$bins
  	
	if (class(multi.hmm)==class.multivariate.hmm) {

  	if (is.null(combinations)) {
  		comb.levels <- levels(bins$combination)
  	} else {
  		comb.levels <- combinations
  	}
  	t <- table(bins$combination) / length(bins)
  	t <- t[names(t) %in% comb.levels]
  	return(t)
  	
	} else if (class(multi.hmm)==class.combined.multivariate.hmm) {
	  
  	if (is.null(combinations)) {
  		comb.levels <- unique(as.vector(sapply(mcols(bins), levels)))
  	} else {
  		comb.levels <- combinations
  	}
    t <- sapply(mcols(bins), function(x) { table(x) / length(bins) })
    t <- t[rownames(t) %in% comb.levels,]
    return(t)
	  
	}
}


#' Transition frequencies of combinatorial states
#'
#' Get a table of transition frequencies between combinatorial states of different \code{\link{multiHMM}}s.
#'
#' @param multi.hmms A list with \code{\link{multiHMM}} objects or a vector with filenames that contain such objects.
#' @param combined.hmm A \code{\link{combinedMultiHMM}} object. If specified, \code{multi.hmms} is ignored.
#' @param zero.states The string(s) which identifies the zero.states.
#' @param combstates Alternative input instead of \code{multi.hmms}: A list of combinatorial state vectors instead of HMMs. If this is specified, \code{multi.hmms} and \code{combined.hmm} will be ignored.
#' @return A data.frame with transition frequencies.
#' @author Aaron Taudt
#' @export
#' @examples 
#'#=== Step 1: Preparation ===
#'## Prepare the file paths. Exchange this with your input and output directories.
#'inputfolder <- system.file("extdata","euratrans", package="chromstaRData")
#'outputfolder <- file.path(tempdir(), 'SHR-BN-example')

#'## Define experiment structure
#'data(experiment_table)
#'print(experiment_table)

#'## Define assembly
#'# This is only necessary if you have BED files, BAM files are handled automatically.
#'# For common assemblies you can also specify them as 'hg19' for example.
#'data(rn4_chrominfo)
#'head(rn4_chrominfo)

#'#=== Step 2: Run Chromstar ===
#'## Run ChromstaR
#'Chromstar(inputfolder, experiment.table=experiment_table,
#'          outputfolder=outputfolder, numCPU=2, binsize=1000, assembly=rn4_chrominfo,
#'          prefit.on.chr='chr12', mode='mark')
#'## Results are stored in 'outputfolder' and can be loaded for further processing
#'list.files(outputfolder)
#'model <- get(load(file.path(outputfolder,'multivariate-combined', 'combined_mode-mark.RData')))

#'#=== Step 3: Analysis ===
#'# Get frequencies
#'transitionFrequencies(combined.hmm=model)
#'
transitionFrequencies <- function(multi.hmms=NULL, combined.hmm=NULL, zero.states="[]", combstates=NULL) {

	if (is.null(combstates)) {
	  if (is.null(combined.hmm)) {
  		## Get combinatorial states in loop to save memory
  		combstates <- list()
  		for (imodel in 1:length(multi.hmms)) {
  			multi.hmm <- suppressMessages( loadHmmsFromFiles(multi.hmms[[imodel]], check.class=class.multivariate.hmm)[[1]] )
  			combstates[[imodel]] <- multi.hmm$bins$combination
  		}
  		names(combstates) <- names(multi.hmms)
	  } else {
	    combined.hmm <- suppressMessages( loadHmmsFromFiles(combined.hmm, check.class=class.combined.multivariate.hmm)[[1]] )
	    combstates <- as.list(mcols(combined.hmm$bins))
	  }
	}
	num.models <- length(combstates)

	### Get transitions for whole genome ###
	combstates$sep <- '<>'
	gentrans <- do.call(paste, combstates)
	freqtrans <- as.data.frame(table(gentrans) / length(gentrans))
	names(freqtrans) <- c('transition','frequency')
	freqtrans <- freqtrans[order(freqtrans$frequency, decreasing=TRUE),]
	# Cumulative frequencies
	freqtrans$cumulative.frequency <- cumsum(freqtrans$frequency)

	### Assigning groups for frequency table ###
	freqtrans$group <- 'other'
	# Stage-specific and constant states
	combstates$sep <- NULL
	levels.combstates <- unique(unlist(lapply(combstates, levels)))
	levels.combstates <- setdiff(levels.combstates, zero.states)
	levels.combstates <- gsub('\\+','\\\\+',levels.combstates)
	levels.combstates <- gsub('\\[','\\\\[', levels.combstates)
	levels.combstates <- gsub('\\]','\\\\]', levels.combstates)
	for (combination in levels.combstates) {
		mask <- intersect(grep(combination, freqtrans$transition), grep(paste(setdiff(levels.combstates,combination), collapse='|'), freqtrans$transition, invert=TRUE))
		freqtrans$group[mask] <- paste0('stage-specific ',gsub('\\\\','',combination))
		mask <- sapply(gregexpr(combination,freqtrans$transition), length) == num.models
		freqtrans$group[mask] <- paste0('constant ', gsub('\\\\','',combination))
	}
	# Zero transitions
	freqtrans.split <- strsplit(sub('<>$','<><>',as.character(freqtrans$transition)),'<>')
	freqtrans.split <- do.call(rbind, freqtrans.split)
	df <- as.data.frame(apply(freqtrans.split, 2, function(x) { x %in% zero.states } ))
	if (ncol(df)==1) {
		iszero <- Reduce('&', df[,1])
	} else {
		iszero <- Reduce('&', as.list(df))
	}
	freqtrans$group[iszero] <- 'zero transition'

	return(freqtrans)
	
	# ### Assigning groups over whole genome ###
	# mapping <- freqtrans$group
	# names(mapping) <- freqtrans$transition
	# gengroups <- mapping[gentrans]
	# gentrans <- data.frame(transition=gentrans, group=gengroups)
	# 
	# 
	# ## Return value ##
	# return(list(frequencies=freqtrans, transitions=gentrans))

}

