unis2pseudomulti <- function(uni.hmm.list) {

	# Load models
	uni.hmm.list <- loadHmmsFromFiles(uni.hmm.list)

	# Extract coordinates and other stuff
	nummod = length(uni.hmm.list)
	bins <- uni.hmm.list[[1]]$bins
	bins$reads <- NULL
	bins$state <- NULL
	numbins = length(uni.hmm.list[[1]]$bins)
	IDs <- unlist(lapply(uni.hmm.list, "[[", "ID"))
	distributions = lapply(uni.hmm.list,"[[","distributions")
	weights = lapply(uni.hmm.list,"[[","weights")

	# Extract the reads
	cat("Extracting reads from uni.hmm.list...")
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	colnames(reads) <- IDs
	for (imod in 1:nummod) {
		reads[,imod] = uni.hmm.list[[imod]]$bins$reads
	}
	maxreads = max(reads)
	bins$reads <- reads
	cat(" done\n")

	## Get combinatorial states
	cat("Getting combinatorial states...")
	combstates.per.bin = combinatorialStates(uni.hmm.list)
	comb.states.table = table(combstates.per.bin)
	comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	numstates <- length(comb.states)
	bins$state <- factor(combstates.per.bin, levels=comb.states)
	cat(" done\n")
	
	## Calculate transition matrix
	cat("Estimating transition matrix...")
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
	cat(" done\n")

	## Return multi.hmm
	multi.hmm <- list()
	multi.hmm$IDs.univariate <- IDs
	multi.hmm$bins <- bins
	## Segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- multi.hmm$bins
		red.gr.list <- GRangesList()
		for (state in comb.states) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		multi.hmm$segments <- red.gr
		seqlengths(multi.hmm$segments) <- seqlengths(multi.hmm$bins)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
	## Parameters
		# Weights
		tstates <- table(combstates.per.bin)
		multi.hmm$weights <- sort(tstates/sum(tstates), decreasing=T)
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
