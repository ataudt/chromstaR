#' Fit a Hidden Markov Model to a ChIP-seq sample.
#'
#' Fit a HMM to a ChIP-seq sample to determine the modification state of genomic regions, e.g. call peaks in the sample.
#'
#' This function is similar to \code{\link{callPeaksUnivariateAllChr}} but allows to pre-fit on a single chromosome instead of the whole genome. This gives a significant performance increase and can help to converge into a better fit in case of unsteady quality for some chromosomes.
#'
#' @param binned.data A \code{\link{GRanges}} object with binned read counts or a file that contains such an object.
#' @param input.data Input control for the experiment. A \code{\link{GRanges}} object with binned read counts or a file that contains such an object.
#' @param prefit.on.chr A chromosome that is used to pre-fit the Hidden Markov Model. Set to \code{NULL} if you don't want to prefit but use the whole genome instead.
#' @param short If \code{TRUE}, the second fitting step is only done with one iteration.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param init One of the following initialization procedures:
#' \describe{
#' \item{\code{standard}}{The negative binomial of state 'unmodified' will be initialized with \code{mean=mean(counts)}, \code{var=var(counts)} and the negative binomial of state 'modified' with \code{mean=mean(counts)+1}, \code{var=var(counts)}. This procedure usually gives the fastest convergence.}
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
#' @param max.mean If \code{mean(counts)>max.mean}, bins with low read counts will be set to 0. This is a workaround to obtain good fits in the case of large bin sizes.
#' @param post.cutoff False discovery rate. code{NULL} means that the state with maximum posterior probability will be chosen, irrespective of its absolute probability (default=code{NULL}).
#' @param control If set to \code{TRUE}, the binned data will be treated as control experiment. That means only state 'zero-inflation' and 'unmodified' will be used in the HMM.
#' @param keep.posteriors If set to \code{TRUE} (default=\code{FALSE}), posteriors will be available in the output. This is useful to change the post.cutoff later, but increases the necessary disk space to store the result.
#' @param keep.densities If set to \code{TRUE} (default=\code{FALSE}), densities will be available in the output. This should only be needed debugging.
#' @param verbosity Verbosity level for the fitting procedure. 0 - No output, 1 - Iterations are printed.
#' @return A \code{\link{uniHMM}} object.
#' @author Aaron Taudt, Maria Colome Tatche
#' @seealso \code{\link{uniHMM}}, \code{\link{callPeaksMultivariate}}
#' @export
#' @examples
#'## Get an example BAM file with ChIP-seq reads
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'data(experiment_table)
#'binned <- binReads(file, experiment.table=experiment_table,
#'                   assembly=rn4_chrominfo, binsizes=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
#'## Check if the fit is ok
#'plotHistogram(hmm)
#'
callPeaksUnivariate <- function(binned.data, input.data=NULL, prefit.on.chr=NULL, short=TRUE, eps=0.1, init="standard", max.time=NULL, max.iter=5000, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff=TRUE, read.cutoff.quantile=1, read.cutoff.absolute=500, max.mean=Inf, post.cutoff=0.5, control=FALSE, keep.posteriors=FALSE, keep.densities=FALSE, verbosity=1) {

    if (class(binned.data) == 'character') { 
        message("Loading file ",binned.data)
        binned.data <- loadHmmsFromFiles(binned.data)[[1]]
    }
    if (!is.null(input.data)) {
        if (class(input.data) == 'character') { 
            message("Loading input file(s) ", paste0(input.data, collapse=', '))
            input.datas <- loadHmmsFromFiles(input.data)
            input.data <- input.datas[[1]]
            if (length(input.datas) > 1) {
                for (i1 in 2:length(input.datas)) {
                    input.data$counts <- input.data$counts + input.datas[[i1]]$counts
                }
            }
        }
    }
    if (!is.null(prefit.on.chr)) {
        if (!prefit.on.chr %in% seqlevels(binned.data)) {
            stop("Could not find chromosome ", prefit.on.chr, " for option 'prefit.on.chr'.")
        }
    }

    if (is.null(prefit.on.chr)) {
        model <- callPeaksUnivariateAllChr(binned.data=binned.data, input.data=input.data, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff=read.cutoff, read.cutoff.quantile=read.cutoff.quantile, read.cutoff.absolute=read.cutoff.absolute, max.mean=max.mean, post.cutoff=post.cutoff, control=control, keep.posteriors=keep.posteriors, keep.densities=FALSE, verbosity=verbosity)
    } else {

        message("Fitting on chromosome ", prefit.on.chr)
        pre.binned.data <- binned.data[seqnames(binned.data)==prefit.on.chr]
        if (!is.null(input.data)) {
            pre.input.data <- input.data[seqnames(input.data)==prefit.on.chr]
        } else {
            pre.input.data <- NULL
        }
        pre.model <- callPeaksUnivariateAllChr(binned.data=pre.binned.data, input.data=pre.input.data, eps=eps, init=init, max.time=max.time, max.iter=max.iter, num.trials=num.trials, eps.try=eps.try, num.threads=num.threads, read.cutoff=read.cutoff, read.cutoff.quantile=read.cutoff.quantile, read.cutoff.absolute=read.cutoff.absolute, max.mean=max.mean, post.cutoff=post.cutoff, control=control, keep.posteriors=FALSE, keep.densities=FALSE, verbosity=verbosity)

        if (short) {
            max.iter=1
        }
        model <- pre.model
        model$bins <- binned.data
        model <- suppressWarnings( callPeaksUnivariateAllChr(binned.data=model, input.data=input.data, eps=eps, max.time=max.time, max.iter=max.iter, num.threads=num.threads, read.cutoff=read.cutoff, read.cutoff.quantile=read.cutoff.quantile, read.cutoff.absolute=read.cutoff.absolute, max.mean=max.mean, post.cutoff=post.cutoff, control=control, keep.posteriors=keep.posteriors, keep.densities=keep.densities, verbosity=verbosity) )

    }

    return(model)

}


#' Fit a Hidden Markov Model to a ChIP-seq sample.
#'
#' Fit a HMM to a ChIP-seq sample to determine the modification state of genomic regions, e.g. call peaks in the sample.
#'
#' The Hidden Markov Model which is used to classify the bins uses 3 states: state 'zero-inflation' with a delta function as emission densitiy (only zero read counts), 'unmodified' and 'modified' with Negative Binomials as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. Please refer to our manuscript at \url{http://dx.doi.org/10.1101/038612} for a detailed description of the method.
#'
#' @param binned.data A \code{\link{GRanges}} object with binned read counts or a file that contains such an object.
#' @param input.data Input control for the experiment. A \code{\link{GRanges}} object with binned read counts or a file that contains such an object.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param init One of the following initialization procedures:
#' \describe{
#' \item{\code{standard}}{The negative binomial of state 'unmodified' will be initialized with \code{mean=mean(counts)}, \code{var=var(counts)} and the negative binomial of state 'modified' with \code{mean=mean(counts)+1}, \code{var=var(counts)}. This procedure usually gives the fastest convergence.}
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
#' @param max.mean If \code{mean(counts)>max.mean}, bins with low read counts will be set to 0. This is a workaround to obtain good fits in the case of large bin sizes.
#' @param post.cutoff False discovery rate. code{NULL} means that the state with maximum posterior probability will be chosen, irrespective of its absolute probability (default=code{NULL}).
#' @param control If set to \code{TRUE}, the binned data will be treated as control experiment. That means only state 'zero-inflation' and 'unmodified' will be used in the HMM.
#' @param keep.posteriors If set to \code{TRUE} (default=\code{FALSE}), posteriors will be available in the output. This is useful to change the post.cutoff later, but increases the necessary disk space to store the result.
#' @param keep.densities If set to \code{TRUE} (default=\code{FALSE}), densities will be available in the output. This should only be needed debugging.
#' @param verbosity Verbosity level for the fitting procedure. 0 - No output, 1 - Iterations are printed.
#' @return A \code{\link{uniHMM}} object.
#' @author Aaron Taudt, Maria Coome Tatche
#' @seealso \code{\link{uniHMM}}, \code{\link{callPeaksMultivariate}}
#' @importFrom stats runif
#' @importFrom S4Vectors Rle runmean
callPeaksUnivariateAllChr <- function(binned.data, input.data=NULL, eps=0.01, init="standard", max.time=NULL, max.iter=NULL, num.trials=1, eps.try=NULL, num.threads=1, read.cutoff=TRUE, read.cutoff.quantile=1, read.cutoff.absolute=500, max.mean=Inf, post.cutoff=0.5, control=FALSE, keep.posteriors=FALSE, keep.densities=FALSE, verbosity=1) {

    ### Define cleanup behaviour ###
    on.exit(.C("C_univariate_cleanup"))

    ### Intercept user input ###
    if (check.positive(eps)!=0) stop("argument 'eps' expects a positive numeric")
    if (is.null(max.time)) { max.time <- -1 } else if (check.nonnegative.integer(max.time)!=0) { stop("argument 'max.time' expects a non-negative integer or NULL") }
    if (is.null(max.iter)) { max.iter <- -1 } else if (check.nonnegative.integer(max.iter)!=0) { stop("argument 'max.iter' expects a non-negative integer or NULL") }
    if (check.positive.integer(num.trials)!=0) stop("argument 'num.trials' expects a positive integer")
    if (!is.null(eps.try)) {
        if (check.positive(eps.try)!=0) stop("argument 'eps.try' expects a positive numeric")
    }
    if (check.positive.integer(num.threads)!=0) stop("argument 'num.threads' expects a positive integer")
    if (post.cutoff>1 | post.cutoff<0) stop("argument 'post.cutoff' has to be between 0 and 1 if specified")
    if (check.logical(keep.posteriors)!=0) stop("argument 'keep.posteriors' expects a logical (TRUE or FALSE)")
    if (check.logical(keep.densities)!=0) stop("argument 'keep.densities' expects a logical (TRUE or FALSE)")
    if (check.integer(verbosity)!=0) stop("argument 'verbosity' expects an integer")
    war <- NULL
    if (is.null(eps.try)) eps.try <- eps
    ## Load binned.data and reuse values if present
    if (class(binned.data) == 'character') { 
        binned.data <- loadHmmsFromFiles(binned.data)[[1]]
    }

    ### Assign variables ###
    if (control) {
        state.labels <- state.labels[1:2] # assigned globally outside this function
        state.distributions <- state.distributions[1:2]
    }
    numstates <- length(state.labels)
    iniproc <- which(init==c("standard","random","empiric")) # transform to int

    ### Assign initial parameters ###
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
        continue.from.univariate.hmm <- TRUE
    } else if (class(binned.data) == 'GRanges') {
        A.initial <- double(length=numstates*numstates)
        proba.initial <- double(length=numstates)
        size.initial <- double(length=numstates)
        prob.initial <- double(length=numstates)
        continue.from.univariate.hmm <- FALSE
    } else {
        stop("argument 'binned.data' expects a GRanges with meta-column 'counts' or a file that contains such an object")
    }

    ## Assign more variables ##
    numbins <- length(binned.data)
    counts <- mcols(binned.data)$counts
    if (keep.densities) { lenDensities <- numbins * numstates } else { lenDensities <- 1 }

    ### Input correction ###
    if (!is.null(input.data)) {
        ptm <- startTimedMessage("Correcting read counts for input ...")
        if (is.character(input.data)) {
            input.data <- loadHmmsFromFiles(input.data)[[1]]
        }
        # Artifacts with super high read count
        index <- which(input.data$counts >= quantile(input.data$counts[input.data$counts>0], 0.999))
        index <- c(index, index-1, index+1) # one neighboring bin to each side
        index <- index[index>0 & index<=length(input.data)] # if we hit chromosome boundaries, bad luck
        counts[index] <- 0
        input.data$counts[index] <- 0
        # Correction factor
        inputf <- S4Vectors::Rle(1, length=length(input.data))
        mean.input.counts <- mean(input.data$counts[input.data$counts>0])
        mask.0 <- input.data$counts > 0
        inputf[mask.0] <- mean.input.counts / as.numeric(S4Vectors::runmean(S4Vectors::Rle(input.data$counts), k=15, endrule='constant'))[mask.0]
        inputf[inputf > 1.5] <- 1
        inputf <- as.numeric(inputf)
        counts <- round(counts * inputf)
        stopTimedMessage(ptm)
    }

    ### Check if there are counts in the data, otherwise HMM will blow up ###
    if (!any(counts!=0)) {
        stop("All counts in data are zero. No univariate HMM done.")
    }

    ### Filter high counts out, makes HMM faster ###
    numfiltered <- 0
    if (continue.from.univariate.hmm) {
        mask <- counts > read.cutoff.absolute
        counts[mask] <- read.cutoff.absolute
        numfiltered <- length(which(mask))
    } else {
        if (read.cutoff) {
            read.cutoff.by.quantile <- as.integer(quantile(counts, read.cutoff.quantile))
            read.cutoff.absolute <- min(read.cutoff.by.quantile, read.cutoff.absolute)
            mask <- counts > read.cutoff.absolute
            counts[mask] <- read.cutoff.absolute
            numfiltered <- length(which(mask))
        }
    }
    if (numfiltered > 0) {
        message("Replaced read counts > ",read.cutoff.absolute, " by ",read.cutoff.absolute," in ",numfiltered," bins to enhance performance (option 'read.cutoff').")
    }

    ### Filter out low read counts that arise when the bin size is larger than optimal (should correct the result to near optimal again) ###
    if (!is.infinite(max.mean)) {
        hist <- graphics::hist(counts[counts>0], breaks=0:max(counts), right=FALSE, plot=FALSE)
        maxhist <- which.max(hist$counts)
        if (maxhist-1 > max.mean) {    # -1 to get from 1-based histogram indices to (0-based) read counts
            # Two empirical rules to remove low counts
            read.counts.to.remove.1 <- which(hist$counts[1:maxhist]<=hist$counts[2]) -1
            minlow <- which.min(hist$counts[2:maxhist])
            read.counts.to.remove <- max(c(read.counts.to.remove.1, 2*minlow))
            index.filtered <- which(counts>0 & counts<=read.counts.to.remove)
            counts[index.filtered] <- 0
            if (length(index.filtered)>0) {
                message(paste0("Replaced read counts <= ",read.counts.to.remove," by 0. This was done because the selected bin size is considered too big for this dataset: The mean of the read counts (zeros removed) is bigger than the specified max.mean = ",max.mean,". Check the fits!"))
            }
        }
    }
    
#     ### Initialization of emission distributions ###
#     reads.greater.0 <- counts[counts>0]
#     mean.counts <- mean(reads.greater.0)
#     var.reads <- var(reads.greater.0)
#     imean <- vector()
#     ivar <- vector()
#     if (init=="standard") {
#         imean['unmodified'] <- mean.counts
#         imean['modified'] <- mean.counts + 1
#         ivar['unmodified'] <- var.reads
#         ivar['modified'] <- var.reads
#         # Make sure variance is bigger than mean for negative binomial
#         for (i1 in 1:length(imean)) {
#             if (imean[i1]>=ivar[i1]) {
#                 ivar[i1] <- imean[i1] + 1
#             }
#         }
#     } else if (init=="random") {
#         for (istate in c('unmodified','modified')) {
#             imean[istate] <- stats::runif(1, min=0, max=10*mean.counts)
#             ivar[istate] <- stats::runif(1, min=imean[istate], max=20*mean.counts)
#         }
#     } else if (init=="empiric") {
#         imean['unmodified'] <- mean.counts/2
#         ivar['unmodified'] <- imean['unmodified']*2
#         imean['modified'] <- mean.counts*2
#         ivar['modified'] <- imean['modified']*2
#     } else if (init %in% seqlevels(binned.data)) {
#     } else if {
#         stop("Unknown initialization procedure: ",init)
#     }


    ## Call univariate in a for loop to enable multiple trials
    if (verbosity==0) {
        ptm <- startTimedMessage("Fitting Hidden Markov Model ...")
    }
    modellist <- list()
    for (i_try in 1:num.trials) {
        if (verbosity>=1) message("------------------------------------ Try ",i_try," of ",num.trials," -------------------------------------")
        hmm <- .C("C_univariate_hmm",
            counts = as.integer(counts), # double* O
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
            size.initial = as.double(size.initial), # double* initial_size
            prob.initial = as.double(prob.initial), # double* initial_prob
            A.initial = as.double(A.initial), # double* initial_A
            proba.initial = as.double(proba.initial), # double* initial_proba
            use.initial.params = as.logical(continue.from.univariate.hmm), # bool* use_initial_params
            num.threads = as.integer(num.threads), # int* num_threads
            error = as.integer(0), # int* error (error handling)
            read.cutoff = as.integer(max(counts)), # int* read_cutoff
            verbosity = as.integer(verbosity) # int* verbosity
        )

        hmm$eps <- eps.try
        if (hmm$loglik.delta > hmm$eps) {
            warning("HMM did not converge in trial run ",i_try,"!")
        }
        # Store model in list
        modellist[[i_try]] <- hmm

        # Set init procedure to random
        iniproc <- which('random'==c("standard","random","empiric")) # transform to int
    }
    if (verbosity==0) {
        stopTimedMessage(ptm)
    }

    if (num.trials > 1) {
        # Select fit with best loglikelihood
        indexmax <- which.max(unlist(lapply(modellist,"[[","loglik")))
        hmm <- modellist[[indexmax]]
        message("Selecting try ", indexmax, " of ", length(modellist), " with best loglikelihood.")
    }

    if (eps != eps.try) {
        if (verbosity==0) {
            ptm <- startTimedMessage("Refining Hidden Markov Model ...")
        }
        # Rerun the HMM with different epsilon and initial parameters from trial run
        if (verbosity>=1) message("------------------------- Rerunning try ",indexmax," with eps = ",eps," -------------------------")
        hmm <- .C("C_univariate_hmm",
            counts = as.integer(counts), # double* O
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
            read.cutoff = as.integer(max(counts)), # int* read_cutoff
            verbosity = as.integer(verbosity) # int* verbosity
        )
        if (verbosity==0) {
            stopTimedMessage(ptm)
        }
    }

    ### Issue warnings ###
    hmm$eps <- eps
    if (hmm$loglik.delta > hmm$eps) {
        war <- warning("HMM did not converge!")
    }
    if (hmm$error == 1) {
        stop("A nan occurred during the Baum-Welch! Parameter estimation terminated prematurely. Check your read counts for very high numbers, they could be the cause for this problem.")
    } else if (hmm$error == 2) {
        stop("An error occurred during the Baum-Welch! Parameter estimation terminated prematurely.")
    }

    ### Make return object ###
        result <- list()
        result$info <- attr(binned.data, 'info')
    ## Get states
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        hmm$posteriors <- matrix(hmm$posteriors, ncol=hmm$num.states)
        colnames(hmm$posteriors) <- state.labels
        threshold <- 1-post.cutoff
        if (control) {
            states <- rep(1, hmm$num.bins)
            states[ hmm$posteriors[,2] >= hmm$posteriors[,1] ] <- 2
        } else {
            states <- rep(NA, hmm$num.bins)
            states[ hmm$posteriors[,3]<threshold & hmm$posteriors[,2]<=hmm$posteriors[,1] ] <- 1
            states[ hmm$posteriors[,3]<threshold & hmm$posteriors[,2]>hmm$posteriors[,1] ] <- 2
            states[ hmm$posteriors[,3]>=threshold ] <- 3
        }
        states <- state.labels[states]
    ## Bin coordinates, posteriors and states
        result$bins <- GRanges(seqnames=seqnames(binned.data),
                                                        ranges=ranges(binned.data),
                                                        counts=hmm$counts,
                                                        state=states) 
        result$bins$posteriors <- hmm$posteriors
        if (!control) {
            result$bins$posterior.modified <- hmm$posteriors[,'modified']
        }
        if (keep.densities) {
            result$bins$densities <- matrix(hmm$densities, ncol=hmm$num.states)
        }
        seqlengths(result$bins) <- seqlengths(binned.data)
        stopTimedMessage(ptm)
    ## Peak score as maximum posterior in that peak
        result$bins$peakScores <- getPeakScore.univariate(result$bins)
    ## Segmentation
        ptm <- startTimedMessage("Making segmentation ...")
        df <- as.data.frame(result$bins)
        red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width',grep('posterior', names(df), value=TRUE), 'counts')))
        result$segments <- methods::as(red.df, 'GRanges')
        seqlengths(result$segments) <- seqlengths(binned.data)[seqlevels(result$segments)]
        if (!keep.posteriors) {
            result$bins$posteriors <- NULL
        }
        stopTimedMessage(ptm)
    ## Peaks
        result$peaks <- result$segments[result$segments$state == 'modified']
        result$peaks$state <- NULL
        result$segments <- NULL
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
        names(result$startProbs) <- state.labels
        result$startProbs.initial <- hmm$proba.initial
        names(result$startProbs.initial) <-state.labels
        # Distributions
        distributions <- data.frame(type=state.distributions, size=hmm$size, prob=hmm$prob, mu=dnbinom.mean(hmm$size,hmm$prob), variance=dnbinom.variance(hmm$size,hmm$prob))
        rownames(distributions) <- state.labels
        result$distributions <- distributions
        distributions.initial <- data.frame(type=state.distributions, size=hmm$size.initial, prob=hmm$prob.initial, mu=dnbinom.mean(hmm$size.initial,hmm$prob.initial), variance=dnbinom.variance(hmm$size.initial,hmm$prob.initial))
        rownames(distributions.initial) <- state.labels
        distributions.initial['zero-inflation',2:5] <- c(0,1,0,0)
        result$distributions.initial <- distributions.initial
        # post.cutoff
        result$post.cutoff <- post.cutoff
    ## Convergence info
        convergenceInfo <- list(eps=eps, loglik=hmm$loglik, loglik.delta=hmm$loglik.delta, num.iterations=hmm$num.iterations, time.sec=hmm$time.sec, max.mean=max.mean, read.cutoff=max(hmm$counts))
        result$convergenceInfo <- convergenceInfo
    ## Add class
        class(result) <- class.univariate.hmm
        
    # Return results
    return(result)
}
