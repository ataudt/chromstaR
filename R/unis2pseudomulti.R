unis2pseudomulti <- function(uni.hmm.list) {

	## Intercept user input
	if (check.univariate.modellist(uni.hmm.list)!=0) {
		cat("Loading univariate HMMs from files ...")
		mlist <- NULL
		for (modelfile in uni.hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		uni.hmm.list <- mlist
		remove(mlist)
		cat(" done\n")
		if (check.univariate.modellist(uni.hmm.list)!=0) stop("argument 'uni.hmm.list' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}

	# Extract coordinates and other stuff
	nummod = length(uni.hmm.list)
	numbins = uni.hmm.list[[1]]$num.bins
	IDs <- unlist(lapply(uni.hmm.list, "[[", "ID"))
	coordinates <- uni.hmm.list[[1]]$coordinates
	seqlengths <- uni.hmm.list[[1]]$seqlengths
	distributions = lapply(uni.hmm.list,"[[","distributions")
	weights = lapply(uni.hmm.list,"[[","weights")

	# Extract the reads
	cat("Extracting reads from uni.hmm.list...")
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	colnames(reads) <- IDs
	for (imod in 1:nummod) {
		reads[,imod] = uni.hmm.list[[imod]]$reads
	}
	maxreads = max(reads)
	cat(" done\n")

	## Get combinatorial states
	cat("Getting combinatorial states...")
	combstates.per.bin = combinatorial.states(uni.hmm.list)
	comb.states.table = table(combstates.per.bin)
	comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	numstates <- length(comb.states)
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
	multi.hmm <- list(coordinates = coordinates,
										IDs.univariate = IDs,
										seqlengths = seqlengths,
										distributions.univariate = distributions,
										weights.univariate = weights,
										comb.states = comb.states,
										states = factor(combstates.per.bin, levels=comb.states),
										reads = reads,
										num.bins = numbins,
										num.states = numstates,
										num.modifications = nummod,
										A.initial = A.estimated)
	class(multi.hmm) <- class.chromstar.multivariate
	return(multi.hmm)

}
