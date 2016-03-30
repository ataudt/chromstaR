#' Combine univariate HMMs to a multivariate HMM
#'
#' Combine multiple \code{\link{chromstaR_univariateHMM}}s to a \code{\link{chromstaR_multivariateHMM}} without running \code{\link{callPeaksMultivariate}}. This should only be done for comparison purposes.
#'
#' Use this function if you want to combine ChIP-seq samples without actually running a multivariate Hidden Markov Model. The resulting object will be of class \code{\link{chromstaR_multivariateHMM}} but will not be truly multivariate.
#'
#' @author Aaron Taudt
#' @param uni.hmm.list A list of \code{\link{chromstaR_univariateHMM}}s, e.g. \code{list(hmm1, hmm2, ...)}.
#' @return A \code{\link{chromstaR_multivariateHMM}} object.
#' @examples
#'## Get example BED-files with ChIP-seq reads for H3K36me3
#' # in 7 different brain tissues (chr22)
#'path.to.example <- system.file(file.path("extdata","brain"), package="chromstaR")
#'bedfiles <- list.files(path.to.example, full=TRUE)
#'## Bin the data into bin size 1000bp and build the univariate HMM
#'uni.HMMs <- list()
#'for (bedfile in bedfiles) {
#'  binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000,
#'                            save.as.RData=FALSE)
#'  uni.HMMs[[bedfile]] <- callPeaksUnivariate(binned.data, ID=basename(bedfile),
#'                                             max.time=30, eps=0.01)
#'}
#'## Combine the univariate HMMs without fitting a multivariate HMM
#'pseudo.multi.HMM <- unis2pseudomulti(uni.HMMs)
#' @export
unis2pseudomulti <- function(uni.hmm.list) {

	# Load models
	uni.hmm.list <- loadHmmsFromFiles(uni.hmm.list)

	# Extract coordinates and other stuff
	nummod = length(uni.hmm.list)
	bins <- uni.hmm.list[[1]]$bins
	bins$counts <- NULL
	bins$state <- NULL
	numbins = length(uni.hmm.list[[1]]$bins)
	IDs <- unlist(lapply(uni.hmm.list, "[[", "ID"))
	names(IDs) <- NULL
	distributions = lapply(uni.hmm.list,"[[","distributions")
	weights = lapply(uni.hmm.list,"[[","weights")

	# Extract the reads
	message("Extracting reads from uni.hmm.list...", appendLF=FALSE)
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	colnames(reads) <- IDs
	for (imod in 1:nummod) {
		reads[,imod] = uni.hmm.list[[imod]]$bins$counts
	}
	maxreads = max(reads)
	bins$counts <- reads
	message(" done")

	## Get combinatorial states
	message("Getting combinatorial states...", appendLF=FALSE)
	combstates.per.bin = combinatorialStates(uni.hmm.list)
	comb.states.table = table(combstates.per.bin)
	comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	numstates <- length(comb.states)
	bins$state <- factor(combstates.per.bin, levels=comb.states)
	message(" done")
	
	## Calculate transition matrix
	message("Estimating transition matrix...", appendLF=FALSE)
	A.estimated = matrix(0, ncol=2^nummod, nrow=2^nummod)
	colnames(A.estimated) = 1:2^nummod-1
	rownames(A.estimated) = 1:2^nummod-1
	for (i1 in 1:(length(combstates.per.bin)-1)) {
		from = combstates.per.bin[i1] + 1
		to = combstates.per.bin[i1+1] + 1
		A.estimated[from,to] = A.estimated[from,to] + 1
	}
	A.estimated = sweep(A.estimated, 1, rowSums(A.estimated), "/")
	# Select only states that are in data
	A.estimated = A.estimated[as.character(comb.states),as.character(comb.states)]
	message(" done")

	## Return multi.hmm
	multi.hmm <- list()
	multi.hmm$IDs <- IDs
	multi.hmm$bins <- bins
	## Segmentation
		message("Making segmentation ...", appendLF=FALSE)
		ptm <- proc.time()
		gr <- multi.hmm$bins
		red.gr.list <- GRangesList()
		for (state in comb.states) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
			if (length(gr)>0) {
				red.gr.list[[length(red.gr.list)+1]] <- red.gr
			}
		}
		red.gr <- GenomicRanges::sort(unlist(red.gr.list))
		multi.hmm$segments <- red.gr
		seqlengths(multi.hmm$segments) <- seqlengths(multi.hmm$bins)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
	## Parameters
		# Weights
		tstates <- table(combstates.per.bin)
		multi.hmm$weights <- sort(tstates/sum(tstates), decreasing=TRUE)
		# Transition matrices
		multi.hmm$transitionProbs <- A.estimated
		# Distributions
		multi.hmm$distributions <- distributions
		names(multi.hmm$distributions) <- IDs
	## Convergence info
		convergenceInfo <- list(eps=Inf, loglik=Inf, loglik.delta=Inf, num.iterations=Inf, time.sec=Inf)
		multi.hmm$convergenceInfo <- convergenceInfo
	## Correlation matrices
# 		multi.hmm$correlation.matrix <- correlationMatrix2use
	## Add class
		class(multi.hmm) <- class.multivariate.hmm

	return(multi.hmm)

}
