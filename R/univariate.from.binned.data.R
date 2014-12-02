univariate.from.binned.data <- function(binned.data, ID, eps=0.001, init="standard", max.time=-1, max.iter=-1, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff.quantile=0.999, max.mean=10, control=FALSE) {

	## Intercept user input
	IDcheck <- ID  #trigger error if not defined
	if (class(binned.data) != 'GRanges') {
		binned.data <- get(load(binned.data))
		if (class(binned.data) != 'GRanges') stop("argument 'binned.data' expects a GRanges with meta-column 'reads' or a file that contains such an object")
	}
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (check.integer(max.time)!=0) stop("argument 'max.time' expects an integer")
	if (check.integer(max.iter)!=0) stop("argument 'max.iter' expects an integer")
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")

	war <- NULL
	if (is.null(eps.try)) eps.try <- eps

	## Assign variables
	if (control) {
		state.labels <- state.labels[1:2] # assigned globally outside this function
	}
	numstates <- length(state.labels)
	numbins <- length(binned.data)
	reads <- mcols(binned.data)$reads
	iniproc <- which(init==c("standard","random","empiric")) # transform to int
	mean.reads <- mean(reads[reads>0])

	# Check if there are reads in the data, otherwise HMM will blow up
	if (!any(reads!=0)) {
		stop("All reads in data are zero. No univariate HMM done.")
	}

	# Filter high reads out, makes HMM faster
	read.cutoff <- quantile(reads, read.cutoff.quantile)
	names.read.cutoff <- names(read.cutoff)
	read.cutoff <- as.integer(read.cutoff)
	mask <- reads > read.cutoff
	reads[mask] <- read.cutoff
	numfiltered <- length(which(mask))
	if (numfiltered > 0) {
		cat(paste0("Replaced read counts > ",read.cutoff," (",names.read.cutoff," quantile) by ",read.cutoff," in ",numfiltered," bins. Set option 'read.cutoff.quantile=1' to disable this filtering. This filtering was done to increase the speed of the HMM and should not affect the results.\n"))
	}

	# Filter out low read counts that arise when the bin size is larger than optimal (should correct the result to near optimal again)
	hist <- hist(reads[reads>0], breaks=0:max(reads), right=FALSE, plot=FALSE)
	maxhist <- which.max(hist$counts)
	if (maxhist-1 > max.mean) {	# -1 to get from 1-based histogram indices to (0-based) read counts
		# Two empirical rules to remove low reads
		read.counts.to.remove.1 <- which(hist$counts[1:maxhist]<=hist$counts[2]) -1
		minlow <- which.min(hist$counts[2:maxhist])
		read.counts.to.remove <- max(c(read.counts.to.remove.1, 2*minlow))
		index.filtered <- which(reads>0 & reads<=read.counts.to.remove)
		reads[index.filtered] <- 0
		if (length(index.filtered)>0) {
			warning(paste0("Replaced read counts <= ",read.counts.to.remove," by 0. This was done because the selected bin size is considered too big for this dataset: The mean of the read counts (zeros removed) is bigger than the specified max.mean = ",max.mean,". However, the results should be good now."))
		}
	}
	
	## Call univariate in a for loop to enable multiple trials
	modellist <- list()
	for (i_try in 1:num.trials) {
		cat("\n\nTry ",i_try," of ",num.trials," ------------------------------\n")
		hmm <- .C("R_univariate_hmm",
			reads = as.integer(reads), # double* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			size = double(length=numstates), # double* size
			prob = double(length=numstates), # double* prob
			num.iterations = as.integer(max.iter), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps.try), # double* eps
			posteriors = double(length=numbins * numstates), # double* posteriors
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			ini.proc = as.integer(iniproc), # int* iniproc
			size.initial = double(length=numstates), # double* initial_size
			prob.initial = double(length=numstates), # double* initial_prob
			A.initial = double(length=numstates*numstates), # double* initial_A
			proba.initial = double(length=numstates), # double* initial_proba
			use.initial.params = as.logical(0), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)

		hmm$eps <- eps.try
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
			reads = as.integer(reads), # double* O
			num.bins = as.integer(numbins), # int* T
			num.states = as.integer(numstates), # int* N
			size = double(length=numstates), # double* size
			prob = double(length=numstates), # double* prob
			num.iterations = as.integer(max.iter), #  int* maxiter
			time.sec = as.integer(max.time), # double* maxtime
			loglik.delta = as.double(eps), # double* eps
			posteriors = double(length=numbins * numstates), # double* posteriors
			A = double(length=numstates*numstates), # double* A
			proba = double(length=numstates), # double* proba
			loglik = double(length=1), # double* loglik
			weights = double(length=numstates), # double* weights
			ini.proc = as.integer(iniproc), # int* iniproc
			size.initial = as.vector(hmm$size), # double* initial_size
			prob.initial = as.vector(hmm$prob), # double* initial_prob
			A.initial = as.vector(hmm$A), # double* initial_A
			proba.initial = as.vector(hmm$proba), # double* initial_proba
			use.initial.params = as.logical(1), # bool* use_initial_params
			num.threads = as.integer(num.threads), # int* num_threads
			error = as.integer(0), # int* error (error handling)
			read.cutoff = as.integer(read.cutoff) # int* read_cutoff
		)
	}

	### Issue warnings ###
	if (num.trials == 1) {
		if (hmm$loglik.delta > eps) {
			war <- warning("HMM did not converge!\n")
		}
	}
	if (hmm$error == 1) {
		stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
	} else if (hmm$error == 2) {
		stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
	}

	### Make return object ###
		result <- list()
		result$ID <- ID
	## Get states
		hmm$posteriors <- matrix(hmm$posteriors, ncol=hmm$num.states)
		colnames(hmm$posteriors) <- paste0("P(",state.labels,")")
		threshold <- 0.5
		states <- rep(NA,hmm$num.bins)
		states[ hmm$posteriors[,3]<=threshold & hmm$posteriors[,2]<=hmm$posteriors[,1] ] <- 1
		states[ hmm$posteriors[,3]<=threshold & hmm$posteriors[,2]>=hmm$posteriors[,1] ] <- 2
		states[ hmm$posteriors[,3]>threshold ] <- 3
		states <- state.labels[states]
	## Bin coordinates and states
		result$bins <- GRanges(seqnames=seqnames(binned.data),
														ranges=ranges(binned.data),
														reads=hmm$reads,
														state=states) 
		seqlengths(result$bins) <- seqlengths(binned.data)
	## Segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- result$bins
		red.gr.list <- GRangesList()
		for (state in state.labels) {
			red.gr <- GenomicRanges::reduce(gr[states==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		result$segments <- red.gr
		seqlengths(result$segments) <- seqlengths(binned.data)
		cat(paste0(" ",round(time[3],2),"s\n"))
	## Parameters
		# Weights
		result$weights <- hmm$weights
		names(result$weights) <- state.labels
		# Transition matrices
		transitionProbs <- matrix(hmm$A, ncol=hmm$num.states, byrow=TRUE)
		rownames(transitionProbs) <- state.labels
		colnames(transitionProbs) <- state.labels
		result$transitionProbs <- transitionProbs
		transitionProbs.initial <- matrix(hmm$A.initial, ncol=hmm$num.states, byrow=TRUE)
		rownames(transitionProbs.initial) <- state.labels
		colnames(transitionProbs.initial) <- state.labels
		result$transitionProbs.initial <- transitionProbs.initial
		# Initial probs
		result$startProbs <- hmm$proba
		names(result$startProbs) <- paste0("P(",state.labels,")")
		result$startProbs.initial <- hmm$proba.initial
		names(result$startProbs.initial) <- paste0("P(",state.labels,")")
		# Distributions
		distributions <- data.frame(type=state.distributions, size=hmm$size, prob=hmm$prob, mu=dnbinom.mean(hmm$size,hmm$prob), variance=dnbinom.variance(hmm$size,hmm$prob))
		rownames(distributions) <- state.labels
		result$distributions <- distributions
		distributions.initial <- data.frame(type=state.distributions, size=hmm$size.initial, prob=hmm$prob.initial, mu=dnbinom.mean(hmm$size.initial,hmm$prob.initial), variance=dnbinom.variance(hmm$size.initial,hmm$prob.initial))
		rownames(distributions.initial) <- state.labels
		distributions.initial['zero-inflation',2:5] <- c(0,1,0,0)
		result$distributions.initial <- distributions.initial
	## Convergence info
		convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec)
		result$convergenceInfo <- convergenceInfo
	## Add class
		class(result) <- class.univariate.hmm
	
	# Return results
	return(result)
}
