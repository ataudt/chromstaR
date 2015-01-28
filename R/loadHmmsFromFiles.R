loadHmmsFromFiles <- function(hmm.list) {

	## Intercept user input
	if (check.univariate.modellist(hmm.list)!=0) {
		message("loading univariate HMMs from files ...", appendLF=F)
		ptm <- proc.time()
		mlist <- NULL
		for (modelfile in hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		hmm.list <- mlist
		remove(mlist)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
		if (check.univariate.modellist(hmm.list)!=0) stop("argument 'hmm.list' expects a list of univariate hmms or a list of files that contain univariate hmms")
	}
	
	return(hmm.list)

}

loadMultiHmmsFromFiles <- function(hmm.list) {

	## Intercept user input
	if (check.multivariate.modellist(hmm.list)!=0) {
		message("loading multivariate HMMs from files ...", appendLF=F)
		ptm <- proc.time()
		mlist <- NULL
		for (modelfile in hmm.list) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		hmm.list <- mlist
		remove(mlist)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
		if (check.multivariate.modellist(hmm.list)!=0) stop("argument 'hmm.list' expects a list of multivariate hmms or a list of files that contain multivariate hmms")
	}
	
	return(hmm.list)

}

