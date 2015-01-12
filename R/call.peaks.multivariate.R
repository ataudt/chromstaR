call.peaks.multivariate <- function(modellist, use.states=NULL, num.states=NULL, eps=0.001, num.threads=1, max.time=-1, max.iter=-1, checkpoint.after.iter=-1, checkpoint.after.time=-1, checkpoint.file=NULL, checkpoint.overwrite=TRUE, checkpoint.use.existing=FALSE, A.initial=NULL) {

	## Intercept user input
	if (check.univariate.modellist(modellist)!=0) {
		cat("Loading univariate HMMs from files ...")
		ptm <- proc.time()
		mlist <- NULL
		for (modelfile in modellist) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		modellist <- mlist
		remove(mlist)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		if (check.univariate.modellist(modellist)!=0) stop("argument 'modellist' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}
	if (!is.null(use.states)) {
		if (check.nonnegative.integer.vector(use.states)!=0) stop("argument 'comb.states' expects a vector of positive integers")
		num.states <- NULL
	}
	if (!is.null(num.states)) {
		if (check.positive.integer(num.states)!=0) stop("argument 'num.states' expects a positive integer")
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.integer(checkpoint.after.iter)!=0) stop("argument 'checkpoint.after.iter' expects an integer")
	if (check.integer(checkpoint.after.time)!=0) stop("argument 'checkpoint.after.time' expects an integer")
	if (check.logical(checkpoint.overwrite)!=0) stop("argument 'checkpoint.overwrite' expects a logical (TRUE or FALSE)")

	## Prepare the HMM
	params <- prepare.multivariate(modellist, use.states=use.states, num.states=num.states, num.threads=num.threads)
	bins <- params$bins
	reads <- params$reads
	numbins <- params$numbins
	nummod <- params$nummod
	comb.states2use <- params$comb.states
	comb.states.per.bin <- params$comb.states.per.bin
	distributions <- params$distributions
	weights <- params$weights
	correlationMatrix2use <- params$correlationMatrix
	correlationMatrixInverse2use <- params$correlationMatrixInverse
	determinant2use <- params$determinant
	usestateTF <- params$usestateTF
	numstates2use <- params$numstates2use
	IDs <- params$IDs
	# Clean up to reduce memory usage
	remove(modellist)
	remove(params)

	## Starting multivariate HMM
	cat("\nStarting multivariate HMM\n")
	cat("Using the following combinatorial states, covering", mean(comb.states.per.bin %in% comb.states2use)*100, "% of the bins:\n", paste(comb.states2use, collapse=" "),"\n")

	# Prepare input for C function
	rs <- unlist(lapply(distributions,"[",2:3,'size'))
	ps <- unlist(lapply(distributions,"[",2:3,'prob'))
	ws1 <- unlist(lapply(weights,"[",1))
	ws2 <- unlist(lapply(weights,"[",2))
	ws3 <- unlist(lapply(weights,"[",3))
	ws <- ws1 / (ws2+ws1)

	# Load checkpoint file if it exists and if desired
	if (is.null(checkpoint.file)) {
		checkpoint.file.exists <- FALSE
	} else {
		if (file.exists(checkpoint.file)) {
			checkpoint.file.exists <- TRUE
		} else {
			checkpoint.file.exists <- FALSE
		}
	}
	if (checkpoint.file.exists & checkpoint.use.existing) {
		cat("Loading checkpoint file ",checkpoint.file,"\n")
		hmm <- get(load(checkpoint.file))
		A.initial <- hmm$transitionProbs
		proba.initial <- hmm$startProbs
		use.initial <- TRUE
	} else {
		if (is.null(A.initial)) {
			A.initial <- double(length=numstates2use*numstates2use)
			proba.initial <- double(length=numstates2use)
			use.initial <- FALSE
		} else {
			proba.initial <- rep(1/numstates2use, numstates2use)
			use.initial <- TRUE
		}
	}
	if (checkpoint.after.iter < 0) {
		checkpoint.after.iter <- max.iter
	}
	if (checkpoint.after.time < 0) {
		checkpoint.after.time <- max.time
	}
	iteration.total <- 0
	time.total <- 0
	repeat{
		# Determine runtime
		if (max.iter > 0) {
			max.iter.temp <- min(checkpoint.after.iter, max.iter-iteration.total)
		} else {
			max.iter.temp <- checkpoint.after.iter
		}
		if (max.time > 0) {
			max.time.temp <- min(checkpoint.after.time, max.time-time.total)
		} else {
			max.time.temp <- checkpoint.after.time
		}
		# Call the C function
		hmm <- .C("R_multivariate_hmm",
			reads = as.integer(as.vector(reads)), # int* multiO
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates2use), # int* N
			num.modifications = as.integer(nummod), # int* Nmod
			comb.states = as.integer(comb.states2use), # int* comb_states
			size = as.double(rs), # double* size
			prob = as.double(ps), # double* prob
			w = as.double(ws), # double* w
			correlation.matrix.inverse = as.double(correlationMatrixInverse2use), # double* cor_matrix_inv
			determinant = as.double(determinant2use), # double* det
			num.iterations = as.integer(max.iter.temp), # int* maxiter
			time.sec = as.integer(max.time.temp), # double* maxtime
			loglik.delta = as.double(eps), # double* eps
			states = integer(length=numbins), # int* states
			A = double(length=numstates2use*numstates2use), # double* A
			proba = double(length=numstates2use), # double* proba
			loglik = double(length=1), # double* loglik
			A.initial = as.double(A.initial), # double* initial A
			proba.initial = as.double(proba.initial), # double* initial proba
			use.initial.params = as.logical(use.initial), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0) # error handling
			)
			
		### Check convergence ###
		if (hmm$error == 1) {
			stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
		} else if (hmm$error == 2) {
			stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
		}

		### Make return object ###
			result <- list()
			result$IDs.univariate <- IDs
		## Bin coordinates and states
			result$bins <- bins
			result$bins$reads <- reads
			result$bins$state <- factor(hmm$states, levels=hmm$comb.states)
		## Segmentation
			cat("Making segmentation ...")
			ptm <- proc.time()
			gr <- result$bins
			red.gr.list <- GRangesList()
			for (state in hmm$comb.states) {
				red.gr <- GenomicRanges::reduce(gr[gr$state==state])
				mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
				red.gr.list[[length(red.gr.list)+1]] <- red.gr
			}
			red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
			result$segments <- red.gr
			seqlengths(result$segments) <- seqlengths(result$bins)
			time <- proc.time() - ptm
			cat(paste0(" ",round(time[3],2),"s\n"))
		## Parameters
			# Weights
			tstates <- table(hmm$states)
			result$weights <- sort(tstates/sum(tstates), decreasing=T)
			# Transition matrices
			result$transitionProbs <- matrix(hmm$A, ncol=numstates2use, byrow=TRUE)
			colnames(result$transitionProbs) <- comb.states2use
			rownames(result$transitionProbs) <- comb.states2use
			result$transitionProbs.initial <- matrix(hmm$A.initial, ncol=numstates2use, byrow=TRUE)
			colnames(result$transitionProbs.initial) <- comb.states2use
			rownames(result$transitionProbs.initial) <- comb.states2use
			# Initial probs
			result$startProbs <- hmm$proba
			names(result$startProbs) <- paste0("P(",comb.states2use,")")
			result$startProbs.initial <- hmm$proba.initial
			names(result$startProbs.initial) <- paste0("P(",comb.states2use,")")
			# Distributions
			result$distributions <- distributions
			names(result$distributions) <- IDs
		## Convergence info
			convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec)
			result$convergenceInfo <- convergenceInfo
		## Correlation matrices
			result$correlation.matrix <- correlationMatrix2use
		## Add class
			class(result) <- class.multivariate.hmm

		# Adjust parameters for the next round
		A.initial <- result$transitionProbs
		proba.initial <- hmm$proba
		use.initial <- TRUE
		iteration.total <- iteration.total + hmm$num.iterations
		result$convergenceInfo$num.iterations <- iteration.total
		time.total <- time.total + hmm$time.sec
		result$convergenceInfo$time.sec <- time.total
		# Test if terminating condition has been reached
		if (result$convergenceInfo$loglik.delta <= result$convergenceInfo$eps | (time.total >= max.time & max.time > 0) | (iteration.total >= max.iter & max.iter > 0)) break
		# Reduce the hmm for checkpointing
		result[c('bins','segments')] <- NULL
		# Save checkpoint
		hmm.checkpoint <- result
		if (checkpoint.overwrite) {
			cat("Saving checkpoint to file ",checkpoint.file,"\n")
			save(hmm.checkpoint, file=checkpoint.file)
		} else {
			cfile <- paste(checkpoint.file,"_iteration_",iteration.total, sep="")
			cat("Saving checkpoint to file ",cfile,"\n")
			save(hmm.checkpoint, file=cfile)
		}
		cat("\nTotal time: ", time.total)
		cat("\nTotal iterations: ", iteration.total)
		cat("\n\nRestarting HMM\n")
	}
	
	## Check convergence
	if (result$convergenceInfo$loglik.delta > result$convergenceInfo$eps) {
		war <- warning("HMM did not converge!\n")
	}
	return(result)

}
