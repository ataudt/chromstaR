#' Load HMMs from files
#'
#' Load \code{\link{uniHMM}} objects from file into a list.
#'
#' @param hmm.list A list of files that contain \code{\link{uniHMM}} objects.
#' @param strict If any of the loaded objects is not a \code{\link{uniHMM}} object, an error (\code{strict=TRUE}) or a warning (\code{strict=FALSE}) will be generated.
#' @return A \code{list()} containing all loaded \code{\link{uniHMM}} objects.
#' @author Aaron Taudt
#' @examples
#'\donttest{
#'## Not run:
#'files.with.hmms <- list.files(folder.with.hmms, full=TRUE,
#'                              pattern='binsize_1000')
#'uni.HMMs <- loadHmmsFromFiles(files.with.hmms)
#'## End(Not run)
#'}
#' @seealso loadMultiHmmsFromFiles
#' @export
loadHmmsFromFiles <- function(hmm.list, strict=FALSE) {

	if (is(hmm.list, class.univariate.hmm)) {
		return(list(hmm.list))
	} else if (is.character(hmm.list)) {
		ptm <- startTimedMessage("Loading univariate HMMs from files ...")
		mlist <- list()
		for (modelfile in hmm.list) {
			mlist[[modelfile]] <- get(load(modelfile))
			if (!is(mlist[[modelfile]], class.univariate.hmm)) {
				if (strict) {
					stopTimedMessage(ptm)
					stop("File ",modelfile," does not contain an ",class.univariate.hmm," object.")
				} else {
					class(mlist[[modelfile]]) <- class.univariate.hmm
					warning("File ",modelfile," does not contain a ",class.univariate.hmm," object. Class attribute corrected.")
				}
			}
		}
		stopTimedMessage(ptm)
		return(mlist)
	} else if (check.univariate.modellist(hmm.list)==0) {
		index <- which(unlist(lapply(hmm.list, function(hmm) { !is(hmm, class.univariate.hmm) })))
		if (length(index)>0) {
			if (strict) {
				stop("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '))
			} else {
				for (ind in index) {
					class(hmm.list[[ind]]) <- class.univariate.hmm
				}
				warning("The following list entries do not contain ",class.univariate.hmm," objects: ", paste(index, collapse=' '),". Class attributes corrected.")
			}
		}
		return(hmm.list)
	} else {
		warning("Returned empty list because 'hmm.list' is neither a character vector of filenames, nor a univariate HMM object, nor a list of univariate HMM objects.")
		return(list())
	}
}

#' Load multivariate HMMs from files
#'
#' Load \code{\link{multiHMM}} objects from file into a list.
#'
#' @param hmm.list A list of files that contain \code{\link{multiHMM}} objects.
#' @param strict If any of the loaded objects is not a \code{\link{multiHMM}} object, an error (\code{strict=TRUE}) or a warning (\code{strict=FALSE}) will be generated.
#' @return A \code{list()} containing all loaded \code{\link{multiHMM}} objects.
#' @author Aaron Taudt
#' @examples
#'\donttest{
#'## Not run:
#'files.with.hmms <- list.files(folder.with.hmms, full=TRUE,
#'                              pattern='binsize_1000')
#'multi.HMMs <- loadMultiHmmsFromFiles(files.with.hmms)
#'## End(Not run)
#'}
#' @seealso loadHmmsFromFiles
#' @export
loadMultiHmmsFromFiles <- function(hmm.list, strict=FALSE) {

	if (is(hmm.list, class.multivariate.hmm)) {
		return(list(hmm.list))
	} else if (is.character(hmm.list)) {
		ptm <- startTimedMessage("Loading multivariate HMMs from files ...")
		mlist <- list()
		for (modelfile in hmm.list) {
			mlist[[modelfile]] <- get(load(modelfile))
			if (!is(mlist[[modelfile]], class.multivariate.hmm)) {
				if (strict) {
					stopTimedMessage(ptm)
					stop("File ",modelfile," does not contain a ",class.multivariate.hmm," object.")
				} else {
					class(modelfile) <- class.multivariate.hmm
					warning("File ",modelfile," does not contain a ",class.multivariate.hmm," object. Class attribute corrected.")
				}
			}
		}
		stopTimedMessage(ptm)
		return(mlist)
	} else if (check.multivariate.modellist(hmm.list)) {
		index <- which(unlist(lapply(hmm.list, function(hmm) { !is(hmm, class.multivariate.hmm) })))
		if (length(index)>0) {
			if (strict) {
				stop("The following list entries do not contain ",class.multivariate.hmm," objects: ", paste(index, collapse=' '))
			} else {
				for (ind in index) {
					class(hmm.list[[ind]]) <- class.multivariate.hmm
				}
				warning("The following list entries do not contain ",class.multivariate.hmm," objects: ", paste(index, collapse=' '),". Class attributes corrected.")
			}
		}
		return(hmm.list)
	} else {
		warning("Returned empty list because 'hmm.list' is neither a character vector of filenames, nor a multivariate HMM object, nor a list of multivariate HMM objects.")
		return(list())
	}
}

