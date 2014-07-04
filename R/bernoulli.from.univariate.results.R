bernoulli.from.univariate.results = function(modellist, use.states=NULL, num.states=NULL, eps=0.001, num.threads=2, max.time=-1, max.iter=-1, output.if.not.converged=FALSE, checkpoint.after.iter=-1, checkpoint.after.time=-1, checkpoint.file="chromStar_checkpoint", checkpoint.reduce=c("coordinates","reads","posteriors","univariate.prob.unmodified"), checkpoint.overwrite=TRUE, checkpoint.use.existing=FALSE, A.initial=NULL) {

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
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")
	if (check.integer(checkpoint.after.iter)!=0) stop("argument 'checkpoint.after.iter' expects an integer")
	if (check.integer(checkpoint.after.time)!=0) stop("argument 'checkpoint.after.time' expects an integer")
	if (check.logical(checkpoint.overwrite)!=0) stop("argument 'checkpoint.overwrite' expects a logical (TRUE or FALSE)")

	## Prepare the HMM
	params = prepare.multivariate(modellist, use.states, num.states)
	coordinates = params$coordinates
	prob.unmodified = params$prob.unmodified
	numbins = params$numbins
	nummod = params$nummod
	comb.states2use = params$comb.states
	distributions = params$distributions
	weights = params$weights
	usestateTF = params$usestateTF
	numstates2use = params$numstates2use
	# Clean up to reduce memory usage
	remove(modellist)

	## Starting multivariate HMM
	cat("\nStarting multivariate HMM\n")
	cat("Using the following states: ",paste(comb.states2use, collapse=" "),"\n")

	# Prepare input for C function

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
		} else {
			proba.initial = rep(1/numstates2use, numstates2use)
			use.initial = TRUE
		}
	}
	if (checkpoint.after.iter < 0) {
		checkpoint.after.iter = max.iter
	}
	if (checkpoint.after.time < 0) {
		checkpoint.after.time = max.time
	}
	iteration.total = 0
	time.total = 0
	repeat{
		# Determine runtime
		if (max.iter > 0) {
			maxiter = min(checkpoint.after.iter, max.iter-iteration.total)
		} else {
			maxiter = checkpoint.after.iter
		}
		if (max.time > 0) {
			maxtime = min(checkpoint.after.time, max.time-time.total)
		} else {
			maxtime = checkpoint.after.time
		}
		# Call the C function
		z = .C("R_multivariate_hmm_productBernoulli",
			as.double(prob.unmodified), # double* multiO
			as.integer(numbins), # int* T
			as.integer(numstates2use), # int* N
			as.integer(nummod), # int* Nmod
			as.integer(comb.states2use), # int* states
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
			as.integer(num.threads), # int* num_threads
			as.integer(0)
			)
			
		# Assign the parameters
		model = NULL
		model$coordinates = coordinates
		model$univariate.prob.unmodified = prob.unmodified
		model$posteriors = matrix(z[[9]], ncol=numstates2use)
		colnames(model$posteriors) = paste("P(state",comb.states2use,")", sep="")
		model$loglik = z[[12]]
		model$iteration = z[[6]]
		model$time.in.sec = z[[7]]
		model$delta.loglik = z[[8]]
		model$epsilon = eps
		model$proba = z[[11]]
		model$A = matrix(z[[10]], ncol=numstates2use, byrow=TRUE)
		colnames(model$A) = comb.states2use
		rownames(model$A) = comb.states2use
		model$distributions.univariate = distributions
		model$weights.univariate = weights
		model$proba.initial = z[[14]]
		model$A.initial = matrix(z[[13]], ncol=numstates2use, byrow=TRUE)
		colnames(model$A.initial) = comb.states2use
		rownames(model$A.initial) = comb.states2use
		model$stateorder = comb.states2use
		model$states = model$stateorder[apply(model$posteriors, 1, which.max)]
		model$num.threads = z[[16]]
		# Adjust parameters for the next round
		A.initial = model$A
		proba.initial = model$proba
		use.initial = TRUE
		iteration.total = iteration.total + model$iteration
		model$iteration = iteration.total
		time.total = time.total + model$time.in.sec
		model$time.in.sec = time.total
		# Test if terminating condition has been reached
		if (model$delta.loglik <= model$epsilon | (time.total >= max.time & max.time > 0) | (iteration.total >= max.iter & max.iter > 0)) break
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
