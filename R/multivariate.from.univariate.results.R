multivariate.from.univariate.results <- function(modellist, use.states=NULL, num.states=NULL, eps=0.001, num.threads=1, max.time=-1, max.iter=-1, output.if.not.converged=FALSE, use.total.correlation=FALSE, checkpoint.after.iter=-1, checkpoint.after.time=-1, checkpoint.file="chromStar_checkpoint", checkpoint.reduce=c("coordinates","reads"), checkpoint.overwrite=TRUE, checkpoint.use.existing=FALSE, A.initial=NULL) {

	## Intercept user input
# 	if (check.univariate.modellist(modellist)!=0) stop("argument 'modellist' expects a list of univariate hmms")
	if (check.univariate.modellist(modellist)!=0) {
		cat("Loading univariate HMMs from files ...")
		mlist <- NULL
		for (modelfile in modellist) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		modellist <- mlist
		remove(mlist)
		cat(" done\n")
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
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")
	if (check.logical(use.total.correlation)!=0) stop("argument 'use.total.correlation' expects a logical (TRUE or FALSE)")
	if (check.integer(checkpoint.after.iter)!=0) stop("argument 'checkpoint.after.iter' expects an integer")
	if (check.integer(checkpoint.after.time)!=0) stop("argument 'checkpoint.after.time' expects an integer")
	if (check.logical(checkpoint.overwrite)!=0) stop("argument 'checkpoint.overwrite' expects a logical (TRUE or FALSE)")

	## Prepare the HMM
	params <- prepare.multivariate(modellist, use.states, num.states)
	coordinates <- params$coordinates
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
	if (use.total.correlation) {
		for (i1 in 1:dim(correlationMatrix2use)[3]) {
			correlationMatrix2use[,,i1] <- params$correlationMatrixAll
			correlationMatrixInverse2use[,,i1] <- params$correlationMatrixAllInverse
			determinant2use[i1] <- params$determinantAll
		}
	}
	usestateTF <- params$usestateTF
	numstates2use <- params$numstates2use
	A.estimated <- params$A.estimated
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
	if (file.exists(checkpoint.file) & checkpoint.use.existing) {
		cat("Loading checkpoint file ",checkpoint.file,"\n")
		hmm <- get(load(checkpoint.file))
		A.initial <- hmm$A
		proba.initial <- hmm$proba
		use.initial <- TRUE
	} else {
		if (is.null(A.initial)) {
			A.initial <- double(length=numstates2use*numstates2use)
			proba.initial <- double(length=numstates2use)
			use.initial <- FALSE
		} else if (A.initial=="estimate") {
			A.initial <- A.estimated
			proba.initial <- rep(1/numstates2use, numstates2use)
			use.initial <- TRUE
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
			maxiter <- min(checkpoint.after.iter, max.iter-iteration.total)
		} else {
			maxiter <- checkpoint.after.iter
		}
		if (max.time > 0) {
			maxtime <- min(checkpoint.after.time, max.time-time.total)
		} else {
			maxtime <- checkpoint.after.time
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
			num.iterations = as.integer(maxiter), # int* maxiter
			time.sec = as.integer(maxtime), # double* maxtime
			loglik.delta = as.double(eps), # double* eps
# 			posteriors = double(length=numbins*numstates2use), # double* posteriors
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
			
		# Add useful entries
		class(hmm) <- class.chromstar.multivariate
		hmm$coordinates <- coordinates
		hmm$reads <- reads # reassign because of matrix layout
# 		hmm$posteriors <- matrix(hmm$posteriors, ncol=numstates2use)
# 		colnames(hmm$posteriors) <- paste("P(state.",comb.states2use,")", sep="")
		hmm$eps <- eps
		hmm$A <- matrix(hmm$A, ncol=numstates2use, byrow=TRUE)
		colnames(hmm$A) <- comb.states2use
		rownames(hmm$A) <- comb.states2use
		hmm$IDs.univariate <- IDs
		hmm$distributions.univariate <- distributions
		hmm$weights.univariate <- weights
		hmm$A.initial <- matrix(hmm$A.initial, ncol=numstates2use, byrow=TRUE)
		colnames(hmm$A.initial) <- comb.states2use
		rownames(hmm$A.initial) <- comb.states2use
# 		hmm$states <- hmm$comb.states[apply(hmm$posteriors, 1, which.max)]
		hmm$correlation.matrix <- correlationMatrix2use
		hmm$correlation.matrix.inverse <- correlationMatrixInverse2use

		# Delete redundant entries
		hmm$size <- NULL
		hmm$prob <- NULL
		hmm$w <- NULL
		hmm$use.initial.params <- NULL

		# Adjust parameters for the next round
		A.initial <- hmm$A
		proba.initial <- hmm$proba
		use.initial <- TRUE
		iteration.total <- iteration.total + hmm$num.iterations
		hmm$num.iterations <- iteration.total
		time.total <- time.total + hmm$time.sec
		hmm$time.sec <- time.total
		# Test if terminating condition has been reached
		if (hmm$loglik.delta <= hmm$eps | (time.total >= max.time & max.time > 0) | (iteration.total >= max.iter & max.iter > 0)) break
		# Reduce the hmm for checkpointing
		hmm[checkpoint.reduce] <- NULL
		# Save checkpoint
		if (checkpoint.overwrite) {
			cat("Saving checkpoint to file ",checkpoint.file,"\n")
			save(hmm, file=checkpoint.file)
		} else {
			cfile <- paste(checkpoint.file,"_iteration_",iteration.total, sep="")
			cat("Saving checkpoint to file ",cfile,"\n")
			save(hmm, file=cfile)
		}
		cat("\nTotal time: ", time.total)
		cat("\nTotal iterations: ", iteration.total)
		cat("\n\nRestarting HMM\n")
	}
	
	# Check convergence
	war <- NULL
	if (hmm$loglik.delta > hmm$eps) {
		war <- warning("HMM did not converge!\n")
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
	} else if (hmm$error == 2) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
	}

	if (!is.null(war)) {
		if (output.if.not.converged == TRUE) {
			return(hmm)
		}
	} else {
		return(hmm)
	}

}
