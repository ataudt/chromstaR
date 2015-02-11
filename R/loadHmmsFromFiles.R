loadHmmsFromFiles <- function(filenames) {

	## Intercept user input
	if (check.univariate.modellist(filenames)!=0) {
		message("loading univariate HMMs from files ...", appendLF=F)
		ptm <- proc.time()
		mlist <- list()
		for (modelfile in filenames) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
		if (check.univariate.modellist(mlist)!=0) stop("argument 'filenames' expects a list of univariate hmms or a list of files that contain univariate hmms")
		return(mlist)
	} else {
		return(filenames)
	}
	

}

loadMultiHmmsFromFiles <- function(filenames) {

	## Intercept user input
	if (check.multivariate.modellist(filenames)!=0) {
		message("loading multivariate HMMs from files ...", appendLF=F)
		ptm <- proc.time()
		mlist <- list()
		for (modelfile in filenames) {
			mlist[[length(mlist)+1]] <- get(load(modelfile))
		}
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
		if (check.multivariate.modellist(mlist)!=0) stop("argument 'filenames' expects a list of multivariate hmms or a list of files that contain multivariate hmms")
		return(mlist)
	} else {
		return(filenames)
	}

}

