#' Fit a Hidden Markov Model to a ChIP-seq sample.
#'
#' Fit a HMM to a ChIP-seq sample to determine the modification state of genomic regions, e.g. call peaks in the sample.
#'
#' This function is similar to \code{\link{callPeaksUnivariate}} but allows to pre-fit on a single chromosome instead of the whole genome. This gives a significant performance increase and can help to converge into a better fit in case of unsteady quality for some chromosomes.
#'
#' @param binned.data A \link{GRanges} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various downstream functions. Could be the file name of the \code{binned.data} for example.
#' @param prefit.on.chr A chromosome that is used to pre-fit the Hidden Markov Model.
#' @param If \code{TRUE}, the second fitting step is only done with one iteration.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param init One of the following initialization procedures:
#' \describe{
#' \item{\code{standard}}{The negative binomial of state 'unmodified' will be initialized with \code{mean=mean(reads)}, \code{var=var(reads)} and the negative binomial of state 'modified' with \code{mean=mean(reads)+1}, \code{var=var(reads)}. This procedure usually gives the fastest convergence.}
#' \item{\code{random}}{Mean and variance of the negative binomials will be initialized with random values (in certain boundaries, see source code). Try this if the \code{'standard'} procedure fails to produce a good fit.}
#' \item{\code{empiric}}{Yet another way to initialize the Baum-Welch. Try this if the other two methods fail to produce a good fit.}
#' }
#' @param max.time The maximum running time in seconds for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default \code{NULL} is no limit.
#' @param max.iter The maximum number of iterations for the Baum-Welch algorithm. The default \code{NULL} is no limit.
#' @param num.trials The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output.
#' @param eps.try If code num.trials is set to greater than 1, \code{eps.try} is used for the trial runs. If unset, \code{eps} is used.
#' @param num.threads Number of threads to use. Setting this to >1 may give increased performance.
#' @param read.cutoff The default (\code{TRUE}) enables filtering of high read counts. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param read.cutoff.quantile A quantile between 0 and 1. Should be near 1. Read counts above this quantile will be set to the read count specified by this quantile. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. If option \code{read.cutoff.absolute} is also specified, the minimum of the resulting cutoff values will be used. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param read.cutoff.absolute Read counts above this value will be set to the read count specified by this value. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. If option \code{read.cutoff.quantile} is also specified, the minimum of the resulting cutoff values will be used. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param max.mean If \code{mean(reads)>max.mean}, bins with low read counts will be set to 0. This is a workaround to obtain good fits in the case of large bin sizes.
#' @param FDR False discovery rate. code{NULL} means that the state with maximum posterior probability will be chosen, irrespective of its absolute probability (default=code{NULL}).
#' @param control If set to \code{TRUE}, the binned data will be treated as control experiment. That means only state 'zero-inflation' and 'unmodified' will be used in the HMM.
#' @param keep.posteriors If set to \code{TRUE} (default=\code{FALSE}), posteriors will be available in the output. This is useful to change the FDR later, but increases the necessary disk space to store the result.
#' @param keep.densities If set to \code{TRUE} (default=\code{FALSE}), densities will be available in the output. This should only be needed debugging.
#' @param verbosity Verbosity level for the fitting procedure. 0 - No output, 1 - Iterations are printed.
#' @author Aaron Taudt, Maria Coome Tatche
#' @seealso \code{\link{chromstaR_univariateHMM}}, \code{\link{callPeaksMultivariate}}
#' @examples
#'## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue
#'bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz",
#'                       package="chromstaR")
#'## Bin the BED file into bin size 1000bp
#'binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000, save.as.RData=FALSE)
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned.data, ID='example_H3K36me3', max.time=60)
#'## Check if the fit is ok
#'plot(hmm, type='histogram')
#' @export
callPeaksUnivariate2 <- function(binned.data, ID, prefit.on.chr, short=TRUE, eps=0.01, init="standard", max.time=NULL, max.iter=NULL, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff=TRUE, read.cutoff.quantile=1, read.cutoff.absolute=500, max.mean=Inf, FDR=0.5, control=FALSE, keep.posteriors=FALSE, keep.densities=FALSE, verbosity=1) {

	if (class(binned.data) == 'character') { 
		message("Loading file ",binned.data)
		binned.data <- get(load(binned.data))
	}

	pre.binned.data <- binned.data[seqnames(binned.data)==prefit.on.chr]
	pre.model <- callPeaksUnivariate(pre.binned.data, ID=ID, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff=read.cutoff, read.cutoff.quantile=read.cutoff.quantile, read.cutoff.absolute=read.cutoff.absolute, max.mean=max.mean, FDR=FDR, control=control, keep.posteriors=FALSE, keep.densities=FALSE, verbosity=verbosity)

	if (short) {
		max.iter=1
	}
	model <- pre.model
	model$bins <- binned.data
	model <- callPeaksUnivariate(model, ID=ID, eps=eps, max.time=max.time, max.iter=max.iter, num.threads=num.threads, read.cutoff=read.cutoff, read.cutoff.quantile=read.cutoff.quantile, read.cutoff.absolute=read.cutoff.absolute, max.mean=max.mean, FDR=FDR, control=control, keep.posteriors=keep.posteriors, keep.densities=keep.densities, verbosity=verbosity)

	return(model)

}


#' Fit a Hidden Markov Model to a ChIP-seq sample.
#'
#' Fit a HMM to a ChIP-seq sample to determine the modification state of genomic regions, e.g. call peaks in the sample.
#'
#' The Hidden Markov Model which is used to classify the bins uses 3 states: state 'zero-inflation' with a delta function as emission densitiy (only zero read counts), 'unmodified' and 'modified' with Negative Binomials as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper TODO:insert paper for a detailed description of the method.
#'
#' @param binned.data A \link{GRanges} object with binned read counts.
#' @param ID An identifier that will be used to identify this sample in various downstream functions. Could be the file name of the \code{binned.data} for example.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param init One of the following initialization procedures:
#' \describe{
#' \item{\code{standard}}{The negative binomial of state 'unmodified' will be initialized with \code{mean=mean(reads)}, \code{var=var(reads)} and the negative binomial of state 'modified' with \code{mean=mean(reads)+1}, \code{var=var(reads)}. This procedure usually gives the fastest convergence.}
#' \item{\code{random}}{Mean and variance of the negative binomials will be initialized with random values (in certain boundaries, see source code). Try this if the \code{'standard'} procedure fails to produce a good fit.}
#' \item{\code{empiric}}{Yet another way to initialize the Baum-Welch. Try this if the other two methods fail to produce a good fit.}
#' }
#' @param max.time The maximum running time in seconds for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default \code{NULL} is no limit.
#' @param max.iter The maximum number of iterations for the Baum-Welch algorithm. The default \code{NULL} is no limit.
#' @param num.trials The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output.
#' @param eps.try If code num.trials is set to greater than 1, \code{eps.try} is used for the trial runs. If unset, \code{eps} is used.
#' @param num.threads Number of threads to use. Setting this to >1 may give increased performance.
#' @param read.cutoff The default (\code{TRUE}) enables filtering of high read counts. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param read.cutoff.quantile A quantile between 0 and 1. Should be near 1. Read counts above this quantile will be set to the read count specified by this quantile. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. If option \code{read.cutoff.absolute} is also specified, the minimum of the resulting cutoff values will be used. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param read.cutoff.absolute Read counts above this value will be set to the read count specified by this value. Filtering very high read counts increases the performance of the Baum-Welch fitting procedure. However, if your data contains very few peaks they might be filtered out. If option \code{read.cutoff.quantile} is also specified, the minimum of the resulting cutoff values will be used. Set \code{read.cutoff=FALSE} to disable this filtering.
#' @param max.mean If \code{mean(reads)>max.mean}, bins with low read counts will be set to 0. This is a workaround to obtain good fits in the case of large bin sizes.
#' @param FDR False discovery rate. code{NULL} means that the state with maximum posterior probability will be chosen, irrespective of its absolute probability (default=code{NULL}).
#' @param control If set to \code{TRUE}, the binned data will be treated as control experiment. That means only state 'zero-inflation' and 'unmodified' will be used in the HMM.
#' @param keep.posteriors If set to \code{TRUE} (default=\code{FALSE}), posteriors will be available in the output. This is useful to change the FDR later, but increases the necessary disk space to store the result.
#' @param keep.densities If set to \code{TRUE} (default=\code{FALSE}), densities will be available in the output. This should only be needed debugging.
#' @param checkpoint.after.iter Write a checkpoint file every n iterations. The default \code{NULL} means no checkpointing for iterations.
#' @param checkpoint.after.time Write a checkpoint file every t seconds. The default \code{NULL} means no checkpointing for time.
#' @param checkpoint.file The name of the checkpoint file that will be written.
#' @param checkpoint.overwrite If set to \code{TRUE}, only one checkpoint file will be written. If set to \code{FALSE}, a new checkpoint file will be written at each checkpoint with the total number of iterations appended.
#' @param checkpoint.use.existing If set to \code{TRUE}, the Baum-Welch fitting procedure will be continued from the HMM in the \code{checkpoint.file}.
#' @param verbosity Verbosity level for the fitting procedure. 0 - No output, 1 - Iterations are printed.
#' @author Aaron Taudt, Maria Coome Tatche
#' @seealso \code{\link{chromstaR_univariateHMM}}, \code{\link{callPeaksMultivariate}}
#' @examples
#'## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue
#'bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz",
#'                       package="chromstaR")
#'## Bin the BED file into bin size 1000bp
#'binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000, save.as.RData=FALSE)
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned.data, ID='example_H3K36me3', max.time=60)
#'## Check if the fit is ok
#'plot(hmm, type='histogram')
#' @export
callPeaksUnivariate <- function(binned.data, ID, eps=0.01, init="standard", max.time=NULL, max.iter=NULL, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff=TRUE, read.cutoff.quantile=1, read.cutoff.absolute=500, max.mean=Inf, FDR=0.5, control=FALSE, keep.posteriors=FALSE, keep.densities=FALSE, checkpoint.after.iter=NULL, checkpoint.after.time=NULL, checkpoint.file=paste0('chromstaR_checkpoint_',ID,'.cpt'), checkpoint.overwrite=TRUE, checkpoint.use.existing=FALSE, verbosity=1) {

	### Define cleanup behaviour ###
	on.exit(.C("R_univariate_cleanup"))

	### Intercept user input ###
	IDcheck <- ID  #trigger error if not defined
	if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
	if (is.null(max.time)) { max.time <- -1 } else if (check.nonnegative.integer(max.time)!=0) { stop("argument 'max.time' expects a non-negative integer") }
	if (is.null(max.iter)) { max.iter <- -1 } else if (check.nonnegative.integer(max.iter)!=0) { stop("argument 'max.iter' expects a non-negative integer") }
	if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
	if (!is.null(eps.try)) {
		if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
	}
	if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
	if (is.null(checkpoint.after.time)) { checkpoint.after.time <- -1 } else if (check.positive.integer(checkpoint.after.time)!=0) { stop("argument 'checkpoint.after.time' expects a positive integer") }
	if (is.null(checkpoint.after.iter)) { checkpoint.after.iter <- -1 } else if (check.positive.integer(checkpoint.after.iter)!=0) { stop("argument 'checkpoint.after.iter' expects a positive integer") }
	if (check.logical(checkpoint.overwrite)!=0) stop("argument 'checkpoint.overwrite' expects a logical (TRUE or FALSE)")
	if (FDR>1 | FDR<0) stop("argument 'FDR' has to be between 0 and 1 if specified")
	if (check.logical(keep.posteriors)!=0) stop("argument 'keep.posteriors' expects a logical (TRUE or FALSE)")
	if (check.logical(keep.densities)!=0) stop("argument 'keep.densities' expects a logical (TRUE or FALSE)")
	if (check.integer(verbosity)!=0) stop("argument 'verbosity' expects an integer")
	war <- NULL
	if (is.null(eps.try)) eps.try <- eps
	## Load binned.data and reuse values if present
	continue.from.univariate.hmm <- FALSE
	if (class(binned.data) == 'character') { 
		message("Loading file ",binned.data)
		binned.data <- get(load(binned.data))
	}
	if (class(binned.data) == class.univariate.hmm) {
		message("Using parameters from univariate HMM")
		hmm <- binned.data
		binned.data <- hmm$bins
		binned.data$state <- NULL
		read.cutoff.absolute <- hmm$convergenceInfo$read.cutoff
		max.mean <- hmm$convergenceInfo$max.mean
		A.initial <- hmm$transitionProbs
		proba.initial <- hmm$startProbs
		size.initial <- hmm$distributions$size
		prob.initial <- hmm$distributions$prob
		use.initial <- TRUE
		continue.from.univariate.hmm <- TRUE
	} else if (class(binned.data) == 'GRanges') {
	} else {
		stop("argument 'binned.data' expects a GRanges with meta-column 'reads' or a file that contains such an object")
	}

	### Assign variables ###
	if (control) {
		state.labels <- state.labels[1:2] # assigned globally outside this function
	}
	numstates <- length(state.labels)
	numbins <- length(binned.data)
	reads <- mcols(binned.data)$reads
	iniproc <- which(init==c("standard","random","empiric")) # transform to int
	if (keep.densities) { lenDensities <- numbins * numstates } else { lenDensities <- 1 }

	### Check if there are reads in the data, otherwise HMM will blow up ###
	if (!any(reads!=0)) {
		stop("All reads in data are zero. No univariate HMM done.")
	}

	### Filter high reads out, makes HMM faster ###
	numfiltered <- 0
	if (!continue.from.univariate.hmm) {
		if (read.cutoff) {
			read.cutoff.by.quantile <- quantile(reads, read.cutoff.quantile)
			read.cutoff.by.quantile <- as.integer(read.cutoff.by.quantile)
			read.cutoff.absolute <- min(read.cutoff.by.quantile, read.cutoff.absolute)
			mask <- reads > read.cutoff.absolute
			reads[mask] <- read.cutoff.absolute
			numfiltered <- length(which(mask))
		}
	} else {
		mask <- reads > read.cutoff.absolute
		reads[mask] <- read.cutoff.absolute
		numfiltered <- length(which(mask))
	}
	if (numfiltered > 0) {
		message("Replaced read counts > ",read.cutoff.absolute, " by ",read.cutoff.absolute," in ",numfiltered," bins. Set option 'read.cutoff=FALSE' to disable this filtering. This filtering was done to increase the speed of the HMM.")
	}

	### Filter out low read counts that arise when the bin size is larger than optimal (should correct the result to near optimal again) ###
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
			message(paste0("Replaced read counts <= ",read.counts.to.remove," by 0. This was done because the selected bin size is considered too big for this dataset: The mean of the read counts (zeros removed) is bigger than the specified max.mean = ",max.mean,". Check the fits!"))
		}
	}
	
# 	### Initialization of emission distributions ###
# 	reads.greater.0 <- reads[reads>0]
# 	mean.reads <- mean(reads.greater.0)
# 	var.reads <- var(reads.greater.0)
# 	imean <- vector()
# 	ivar <- vector()
# 	if (init=="standard") {
# 		imean['unmodified'] <- mean.reads
# 		imean['modified'] <- mean.reads + 1
# 		ivar['unmodified'] <- var.reads
# 		ivar['modified'] <- var.reads
# 		# Make sure variance is bigger than mean for negative binomial
# 		for (i1 in 1:length(imean)) {
# 			if (imean[i1]>=ivar[i1]) {
# 				ivar[i1] <- imean[i1] + 1
# 			}
# 		}
# 	} else if (init=="random") {
# 		for (istate in c('unmodified','modified')) {
# 			imean[istate] <- runif(1, min=0, max=10*mean.reads)
# 			ivar[istate] <- runif(1, min=imean[istate], max=20*mean.reads)
# 		}
# 	} else if (init=="empiric") {
# 		imean['unmodified'] <- mean.reads/2
# 		ivar['unmodified'] <- imean['unmodified']*2
# 		imean['modified'] <- mean.reads*2
# 		ivar['modified'] <- imean['modified']*2
# 	} else if (init %in% seqlevels(binned.data)) {
# 	} else if {
# 		stop("Unknown initialization procedure: ",init)
# 	}

	if (num.trials == 1) {
		### Load checkpoint file if it exists and if desired ###
		if (file.exists(checkpoint.file) & checkpoint.use.existing) {
			message("Loading checkpoint file ",checkpoint.file)
			hmm <- get(load(checkpoint.file))
			message("Using parameters from checkpoint file ",checkpoint.file)
			A.initial <- hmm$transitionProbs
			proba.initial <- hmm$startProbs
			size.initial <- hmm$distributions$size
			prob.initial <- hmm$distributions$prob
			use.initial <- TRUE
		} else if (!continue.from.univariate.hmm) {
			A.initial <- double(length=numstates*numstates)
			proba.initial <- double(length=numstates)
			size.initial <- double(length=numstates)
			prob.initial <- double(length=numstates)
			use.initial <- FALSE
		}
		if (checkpoint.after.iter < 0) {
			checkpoint.after.iter <- max.iter
		}
		if (checkpoint.after.time < 0) {
			checkpoint.after.time <- max.time
		}
		### Run univariate HMM ###
		iteration.total <- 0
		time.total <- 0
		if (verbosity>=1) message("------------------------------------ Try ",1," of ",1," -------------------------------------")
		repeat {
			## Determine runtime
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
			## Call the Baum-Welch
			hmm <- .C("R_univariate_hmm",
				reads = as.integer(reads), # double* O
				num.bins = as.integer(numbins), # int* T
				num.states = as.integer(numstates), # int* N
				size = double(length=numstates), # double* size
				prob = double(length=numstates), # double* prob
				num.iterations = as.integer(max.iter.temp), #  int* maxiter
				time.sec = as.integer(max.time.temp), # double* maxtime
				loglik.delta = as.double(eps.try), # double* eps
				posteriors = double(length=numbins * numstates), # double* posteriors
				densities = double(length=lenDensities), # double* densities
				keep.densities = as.logical(keep.densities), # bool* keep_densities
				A = double(length=numstates*numstates), # double* A
				proba = double(length=numstates), # double* proba
				loglik = double(length=1), # double* loglik
				weights = double(length=numstates), # double* weights
				ini.proc = as.integer(iniproc), # int* iniproc
				size.initial = as.double(size.initial), # double* initial_size
				prob.initial = as.double(prob.initial), # double* initial_prob
				A.initial = as.double(A.initial), # double* initial_A
				proba.initial = as.double(proba.initial), # double* initial_proba
				use.initial.params = as.logical(use.initial), # bool* use_initial_params
				num.threads = as.integer(num.threads), # int* num_threads
				error = as.integer(0), # int* error (error handling)
				read.cutoff = as.integer(max(reads)), # int* read_cutoff
				verbosity = as.integer(verbosity) # int* verbosity
			)
			## Adjust parameters for the next round
			A.initial <- hmm$A
			proba.initial <- hmm$proba
			size.initial <- hmm$size
			prob.initial <- hmm$prob
			use.initial <- TRUE
			iteration.total <- iteration.total + hmm$num.iterations
			hmm$num.iterations <- iteration.total
			time.total <- time.total + hmm$time.sec
			hmm$time.sec <- time.total
			if (hmm$loglik.delta <= eps | (time.total >= max.time & max.time >= 0) | (iteration.total >= max.iter & max.iter >= 0)) break

			### Save checkpoint ###
				result <- list()
				result$ID <- ID
			## Parameters
				# Weights
				result$weights <- hmm$weights
				names(result$weights) <- state.labels
				# Transition matrices
				transitionProbs <- matrix(hmm$A, ncol=hmm$num.states)
				rownames(transitionProbs) <- state.labels
				colnames(transitionProbs) <- state.labels
				result$transitionProbs <- transitionProbs
				transitionProbs.initial <- matrix(hmm$A.initial, ncol=hmm$num.states)
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
				# FDR
				result$FDR <- FDR
			## Convergence info
				convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec, read.cutoff=max(hmm$reads))
				result$convergenceInfo <- convergenceInfo
			## Add class
				class(result) <- class.univariate.hmm
			## Save to file
				hmm.checkpoint <- result
				if (checkpoint.overwrite) {
					message("Saving checkpoint to file ",checkpoint.file)
					save(hmm.checkpoint, file=checkpoint.file)
				} else {
					cfile <- paste(checkpoint.file,"_iteration_",iteration.total, sep="")
					message("Saving checkpoint to file ",cfile)
					save(hmm.checkpoint, file=cfile)
				}

			message("Total time: ", time.total)
			message("Total iterations: ", iteration.total)
			message("Restarting HMM")
		}

	} else {

		## Call univariate in a for loop to enable multiple trials
		modellist <- list()
		for (i_try in 1:num.trials) {
			if (verbosity>=1) message("------------------------------------ Try ",i_try," of ",num.trials," -------------------------------------")
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
				densities = double(length=lenDensities), # double* densities
				keep.densities = as.logical(keep.densities), # bool* keep_densities
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
				read.cutoff = as.integer(max(reads)), # int* read_cutoff
				verbosity = as.integer(verbosity) # int* verbosity
			)

			hmm$eps <- eps.try
			if (hmm$loglik.delta > hmm$eps) {
				warning("HMM did not converge in trial run ",i_try,"!\n")
			}
			# Store model in list
			modellist[[i_try]] <- hmm

			# Set init procedure to random
			iniproc <- which('random'==c("standard","random","empiric")) # transform to int
		}

		# Select fit with best loglikelihood
		indexmax <- which.max(unlist(lapply(modellist,"[[","loglik")))
		hmm <- modellist[[indexmax]]

		# Rerun the HMM with different epsilon and initial parameters from trial run
		if (verbosity>=1) message("------------------------- Rerunning try ",indexmax," with eps = ",eps," -------------------------")
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
			densities = double(length=lenDensities), # double* densities
			keep.densities = as.logical(keep.densities), # bool* keep_densities
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
			read.cutoff = as.integer(max(reads)), # int* read_cutoff
			verbosity = as.integer(verbosity) # int* verbosity
		)
	}

	### Issue warnings ###
	hmm$eps <- eps
	if (hmm$loglik.delta > hmm$eps) {
		war <- warning("HMM did not converge!\n")
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
		message("Calculating states from posteriors ...", appendLF=F); ptm <- proc.time()
		hmm$posteriors <- matrix(hmm$posteriors, ncol=hmm$num.states)
		colnames(hmm$posteriors) <- paste0("P(",state.labels,")")
		threshold <- 1-FDR
		states <- rep(NA,hmm$num.bins)
		states[ hmm$posteriors[,3]<threshold & hmm$posteriors[,2]<=hmm$posteriors[,1] ] <- 1
		states[ hmm$posteriors[,3]<threshold & hmm$posteriors[,2]>hmm$posteriors[,1] ] <- 2
		states[ hmm$posteriors[,3]>=threshold ] <- 3
		states <- state.labels[states]
	## Bin coordinates, posteriors and states
		result$bins <- GRanges(seqnames=seqnames(binned.data),
														ranges=ranges(binned.data),
														reads=hmm$reads,
														state=states) 
		result$bins$score <- apply(hmm$posteriors, 1, max)
		result$bins$posteriors <- hmm$posteriors
		if (keep.densities) {
			result$bins$densities <- matrix(hmm$densities, ncol=hmm$num.states)
		}
		seqlengths(result$bins) <- seqlengths(binned.data)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	## Segmentation
		message("Making segmentation ...", appendLF=F); ptm <- proc.time()
		gr <- result$bins
		red.df <- suppressMessages(collapseBins(as.data.frame(gr), column2collapseBy='state', columns2average=c('reads','score','posteriors.P.modified.'), columns2drop=c('width','posteriors.P.zero.inflation.','posteriors.P.unmodified.')))
		red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], mean.reads=red.df[,'mean.reads'], state=red.df[,'state'], score=red.df[,'mean.score'], mean.posterior.modified=red.df[,'mean.posteriors.P.modified.'])
		result$segments <- red.gr
		seqlengths(result$segments) <- seqlengths(binned.data)
		if (!keep.posteriors) {
			result$bins$posteriors <- NULL
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	## Parameters
		# Weights
		result$weights <- hmm$weights
		names(result$weights) <- state.labels
		# Transition matrices
		transitionProbs <- matrix(hmm$A, ncol=hmm$num.states)
		rownames(transitionProbs) <- state.labels
		colnames(transitionProbs) <- state.labels
		result$transitionProbs <- transitionProbs
		transitionProbs.initial <- matrix(hmm$A.initial, ncol=hmm$num.states)
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
		# FDR
		result$FDR <- FDR
	## Convergence info
		convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec, max.mean=max.mean, read.cutoff=max(hmm$reads))
		result$convergenceInfo <- convergenceInfo
	## Add class
		class(result) <- class.univariate.hmm
	
	# Return results
	return(result)
}
