univariate.from.binned.data <- function(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE) {

	## Intercept user input
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.it)!=0) stop("argument 'max.it' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (check.logical(output.if.not.converged)!=0) stop("argument 'output.if.not.converged' expects a logical (TRUE or FALSE)")


	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	names(binned.data) <- binned.data.names # defined globally outside this function

	## Assign variables
# 	state.labels # assigned globally outside this function
	numstates <- length(state.labels)
	numbins <- length(binned.data$reads)

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(binned.data$reads!=0)) {
		stop("All reads in data are zero. No univariate HMM done.")
	}

	# Filter high reads out, makes HMM faster
	if (filter.reads) {
		limit <- 10*ceiling(var(binned.data$reads))
		mask <- binned.data$reads > limit
		binned.data$reads[mask] <- limit
		numfiltered <- length(which(mask))
		if (numfiltered > 0) {
			warning(paste("There are very high read counts (probably artificial) in your data. Replaced read counts > ",limit," (10*variance) by ",limit," in ",numfiltered," bins. Set option 'filter.reads=FALSE' to disable this filtering.", sep=""))
		}
	}
	
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(binned.data$reads), # double* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			r = double(length=numstates), # double* r
			p = double(length=numstates), # double* p
			num.iterations = as.integer(max.it), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps.try), # double* eps
			posteriors = double(length=numbins * numstates), # double* posteriors
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			r.initial = double(length=numstates), # double* initial_r
			p.initial = double(length=numstates), # double* initial_p
			A.initial = double(length=numstates*numstates), # double* initial_A
			proba.initial = double(length=numstates), # double* initial_proba
			use.initial.params = as.logical(0), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0) # error handling
		)

		names(hmm$weights) <- state.labels
		hmm$eps <- eps.try
		hmm$A <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
		rownames(hmm$A) <- state.labels
		colnames(hmm$A) <- state.labels
		hmm$distributions <- cbind(r=hmm$r, p=hmm$p, mean=fmean(hmm$r,hmm$p), variance=fvariance(hmm$r,hmm$p))
		rownames(hmm$distributions) <- state.labels
		hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
		rownames(hmm$A.initial) <- state.labels
		colnames(hmm$A.initial) <- state.labels
		hmm$distributions.initial <- cbind(r=hmm$r.initial, p=hmm$p.initial, mean=fmean(hmm$r.initial,hmm$p.initial), variance=fvariance(hmm$r.initial,hmm$p.initial))
		rownames(hmm$distributions.initial) <- state.labels
		if (num.trials > 1) {
			if (hmm$loglik.delta > hmm$eps) {
				warning("HMM did not converge in trial run ",i_try,"!\n")
			}
			# Store model in list
			hmm$posteriors <- NULL
			hmm$reads <- NULL
			modellist[[i_try]] <- hmm
		}
	}

	if (num.trials > 1) {

		# Select fit with best loglikelihood
		indexmax <- which.max(unlist(lapply(modellist,"[[","loglik")))
		hmm <- modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		cat("\n\nRerunning try ",indexmax," with eps =",eps,"--------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads <- as.integer(binned.data$reads), # double* O
			num.bins <- as.integer(numbins), # int* T
			num.states <- as.integer(numstates), # int* N
			r <- double(length=numstates), # double* r
			p <- double(length=numstates), # double* p
			num.iterations <- as.integer(max.it), #  int* maxiter
			time.sec <- as.integer(max.time), # double* maxtime
			loglik.delta <- as.double(eps), # double* eps
			posteriors <- double(length=numbins * numstates), # double* posteriors
			A <- double(length=numstates*numstates), # double* A
			proba <- double(length=numstates), # double* proba
			loglik <- double(length=1), # double* loglik
			weights <- double(length=numstates), # double* weights
			r.initial <- as.vector(hmm$distributions[,'r']), # double* initial_r
			p.initial <- as.vector(hmm$distributions[,'p']), # double* initial_p
			A.initial <- as.vector(hmm$A), # double* initial_A
			proba.initial <- as.vector(hmm$proba), # double* initial_proba
			use.initial.params <- as.logical(1), # bool* use_initial_params
			num.threads <- as.integer(num.threads), # int* num_threads
			error = as.integer(0) # error handling
		)
	}

	# Add useful entries
	names(hmm$weights) <- state.labels
	hmm$coordinates <- binned.data[,coordinate.names]
	hmm$posteriors <- matrix(hmm$posteriors, ncol=hmm$num.states)
	colnames(hmm$posteriors) <- paste("P(",state.labels,")", sep="")
	class(hmm) <- class.chromstar.univariate
	class(hmm) <- class.chromstar.univariate
	hmm$states <- get.states(hmm)
	hmm$eps <- eps
	hmm$A <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A) <- state.labels
	colnames(hmm$A) <- state.labels
	hmm$distributions <- cbind(r=hmm$r, p=hmm$p, mean=fmean(hmm$r,hmm$p), variance=fvariance(hmm$r,hmm$p))
	rownames(hmm$distributions) <- state.labels
	hmm$A.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
	rownames(hmm$A.initial) <- state.labels
	colnames(hmm$A.initial) <- state.labels
	hmm$distributions.initial <- cbind(r=hmm$r.initial, p=hmm$p.initial, mean=fmean(hmm$r.initial,hmm$p.initial), variance=fvariance(hmm$r.initial,hmm$p.initial))
	rownames(hmm$distributions.initial) <- state.labels

	# Delete redundant entries
	hmm$r <- NULL
	hmm$p <- NULL
	hmm$r.initial <- NULL
	hmm$p.initial <- NULL
	hmm$use.initial.params <- NULL

	# Issue warnings
	if (num.trials == 1) {
		if (hmm$loglik.delta > hmm$eps) {
			war <- warning("HMM did not converge!\n")
		}
	}
	if (hmm$error == -1) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
	}

	# Return results
	if (!is.null(war)) {
		if (output.if.not.converged == TRUE) {
			return(hmm)
		}
	} else {
		return(hmm)
	}
}
