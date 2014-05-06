multivariate.from.univariate.results = function(modellist, use.states=NULL, num.states=NULL, eps=0.001, num.threads=1, max.time=-1, max.it=-1, output.if.not.converged=FALSE, use.total.correlation=FALSE, checkpoint.after.it=-1, checkpoint.after.time=-1, checkpoint.file="chromStar_checkpoint", checkpoint.reduce=c("coordinates","reads","posteriors"), checkpoint.overwrite=TRUE, checkpoint.use.existing=FALSE, A.initial=NULL) {

	## Intercept user input
	if (check.univariate.modellist(modellist)!=0) stop("argument 'modellist' expects a list of univariate models")
	if (!is.null(use.states)) {
		if (check.nonnegative.integer.vector(use.states)!=0) stop("argument 'comb.states' expects a vector of positive integers")
		num.states = NULL
	}
	if (!is.null(num.states)) {
		if (check.positive.integer(num.states)!=0) stop("argument 'num.states' expects a positive integer")
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.it)!=0) stop("argument 'max.it' expects an integer")
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")
	if (check.logical(use.total.correlation)!=0) stop("argument 'use.total.correlation' expects a logical (TRUE or FALSE)")
	if (check.integer(checkpoint.after.it)!=0) stop("argument 'checkpoint.after.it' expects an integer")
	if (check.integer(checkpoint.after.time)!=0) stop("argument 'checkpoint.after.time' expects an integer")
	if (check.logical(checkpoint.overwrite)!=0) stop("argument 'checkpoint.overwrite' expects a logical (TRUE or FALSE)")

	## Prepare the HMM
	params = prepare.multivariate(modellist, use.states, num.states)
	coordinates = params$coordinates
	reads = params$reads
	numbins = params$numbins
	nummod = params$nummod
	comb.states2use = params$comb.states
	distributions = params$distributions
	softweights = params$softweights
	correlationMatrix2use = params$correlationMatrix
	correlationMatrixInverse2use = params$correlationMatrixInverse
	determinant2use = params$determinant
	if (use.total.correlation) {
		for (i1 in 1:dim(correlationMatrix2use)[3]) {
			correlationMatrix2use[,,i1] = params$correlationMatrixAll
			correlationMatrixInverse2use[,,i1] = params$correlationMatrixAllInverse
			determinant2use[i1] = params$determinantAll
		}
	}
	usestateTF = params$usestateTF
	numstates2use = params$numstates2use
	# Clean up to reduce memory usage
	remove(modellist)

	## Starting multivariate HMM
	cat("\nStarting multivariate HMM\n")
	cat("Using the following states: ",paste(comb.states2use, collapse=" "),"\n")

	# Prepare input for C function
	rs = unlist(lapply(distributions,"[",2:3,'r'))
	ps = unlist(lapply(distributions,"[",2:3,'p'))
	ws = unlist(lapply(softweights,"[",1))

	# Load checkpoint file if it exists and if desired
	if (file.exists(checkpoint.file) & checkpoint.use.existing) {
		cat("Loading checkpoint file ",checkpoint.file,"\n")
		model = get(load(checkpoint.file))
		A.initial = model$A
		proba.initial = model$proba
		use.initial = TRUE
	} else {
		if (is.null(A.initial)) {
			A.initial = double(length=numstates2use*numstates2use)
			proba.initial = double(length=numstates2use)
			use.initial = FALSE
		} else if (A.initial=="estimate") {
			A.initial = params$A.estimated
			proba.initial = rep(1/numstates2use, numstates2use)
			use.initial = TRUE
		} else {
			proba.initial = rep(1/numstates2use, numstates2use)
			use.initial = TRUE
		}
	}
	if (checkpoint.after.it < 0) {
		checkpoint.after.it = max.it
	}
	if (checkpoint.after.time < 0) {
		checkpoint.after.time = max.time
	}
	iteration.total = 0
	time.total = 0
	repeat{
		# Determine runtime
		if (max.it > 0) {
			maxiter = min(checkpoint.after.it, max.it-iteration.total)
		} else {
			maxiter = checkpoint.after.it
		}
		if (max.time > 0) {
			maxtime = min(checkpoint.after.time, max.time-time.total)
		} else {
			maxtime = checkpoint.after.time
		}
		# Call the C function
		z = .C("R_multivariate_hmm",
			as.integer(as.vector(reads)), # int* multiO
			as.integer(numbins), # int* T
			as.integer(numstates2use), # int* N
			as.integer(nummod), # int* Nmod
			as.integer(comb.states2use), # int* states
			as.double(rs), # double* r
			as.double(ps), # double* p
			as.double(ws), # double* w
			as.double(correlationMatrixInverse2use), # double* cor_matrix_inv
			as.double(determinant2use), # double* det
			as.integer(maxiter), # int* maxiter
			as.integer(maxtime), # double* maxtime
			as.double(eps), # double* eps
			double(length=numbins*numstates2use), # double* post
			double(length=numstates2use*numstates2use), # double* A
			double(length=numstates2use), # double* proba
			double(length=1), # double* loglik
			as.double(A.initial), # double* initial A
			as.double(proba.initial), # double* initial proba
			as.logical(use.initial), # bool* use_initial_params
			as.integer(num.threads) # int* num_threads
			)
			
		# Assign the parameters
		model = NULL
		model$coordinates = coordinates
		model$reads = reads
		model$posteriors = matrix(z[[14]], ncol=numstates2use)
		colnames(model$posteriors) = paste("P(state",comb.states2use,")", sep="")
		model$loglik = z[[17]]
		model$iteration = z[[11]]
		model$time.in.sec = z[[12]]
		model$delta.loglik = z[[13]]
		model$epsilon = eps
		model$proba = z[[16]]
		model$A = matrix(z[[15]], ncol=numstates2use, byrow=TRUE)
		colnames(model$A) = comb.states2use
		rownames(model$A) = comb.states2use
		model$distributions.univariate = distributions
		model$softweights.univariate = softweights
		model$proba.initial = z[[19]]
		model$A.initial = matrix(z[[18]], ncol=numstates2use, byrow=TRUE)
		colnames(model$A.initial) = comb.states2use
		rownames(model$A.initial) = comb.states2use
		model$stateorder = comb.states2use
		model$states = model$stateorder[apply(model$posteriors, 1, which.max)]
		model$correlationMatrix = correlationMatrix2use
		model$correlationMatrixInverse = correlationMatrixInverse2use
		model$determinant = determinant2use
		model$num.threads = z[[21]]
		# Adjust parameters for the next round
		A.initial = model$A
		proba.initial = model$proba
		use.initial = TRUE
		iteration.total = iteration.total + model$iteration
		model$iteration = iteration.total
		time.total = time.total + model$time.in.sec
		model$time.in.sec = time.total
		# Test if terminating condition has been reached
		if (model$delta.loglik <= model$epsilon | (time.total >= max.time & max.time > 0) | (iteration.total >= max.it & max.it > 0)) break
		# Reduce the model for checkpointing
		model[checkpoint.reduce] = NULL
		# Save checkpoint
		if (checkpoint.overwrite) {
			cat("Saving checkpoint to file ",checkpoint.file,"\n")
			save(model, file=checkpoint.file)
		} else {
			cfile = paste(checkpoint.file,"_iteration_",iteration.total, sep="")
			cat("Saving checkpoint to file ",cfile,"\n")
			save(model, file=cfile)
		}
		cat("\nTotal time: ", time.total)
		cat("\nTotal iterations: ", iteration.total)
		cat("\n\nRestarting HMM\n")
	}
	
	# Check convergence
	war = NULL
	if (model$delta.loglik > model$epsilon) {
		war = warning("HMM did not converge!\n")
	}
	class(model) = "chromStar.multivariate.model"

	if (!is.null(war)) {
		if (output.if.not.converged == TRUE) {
			return(model)
		}
	} else {
		return(model)
	}

}
