loadHmmsFromFiles <- function(hmm.list) {

	## Intercept user input
	if (check.univariate.modellist(hmm.list)!=0) {
		cat("loading univariate HMMs from files ...")
		ptm <- proc.time()
		mlist <- NULL
		for (modelfile in hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		hmm.list <- mlist
		remove(mlist)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		if (check.univariate.modellist(hmm.list)!=0) stop("argument 'hmm.list' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}
	
	return(hmm.list)

}

loadMultiHmmsFromFiles <- function(hmm.list) {

	## Intercept user input
	if (check.multivariate.modellist(hmm.list)!=0) {
		cat("loading multivariate HMMs from files ...")
		ptm <- proc.time()
		mlist <- NULL
		for (modelfile in hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		hmm.list <- mlist
		remove(mlist)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		if (check.multivariate.modellist(hmm.list)!=0) stop("argument 'hmm.list' expects a list of multivariate hmms or a list of files that contain multivariate hmms")
	}
	
	return(hmm.list)

}

