#' Merge several \code{\link{multiHMM}}s into one object
#'
#' Merge several \code{\link{multiHMM}}s into one object. This can be done to merge fits for separate chromosomes into one object for easier handling. Merging will only be done if all models have the same IDs.
#'
#' @author Aaron Taudt
#' @param multi.hmm.list A list of \code{\link{multiHMM}} objects or a character vector of files that contain such objects.
#' @param filename The file name where the merged object will be stored. If \code{filename} is not specified, a \code{\link{multiHMM}} is returned.
#' @return A \code{\link{multiHMM}} object or NULL, depending on option \code{filename}.
mergeChroms <- function(multi.hmm.list, filename=NULL) {

    ## Check user input
    multi.hmm.list <- loadHmmsFromFiles(multi.hmm.list, check.class=class.multivariate.hmm)
    ## Check if all models have the same ID
    same.IDs <- Reduce('|', unlist(lapply(multi.hmm.list, function(x) { x$info$ID == multi.hmm.list[[1]]$info$ID })))
    if (!is.null(same.IDs)) {
        if (!same.IDs) {
            stop("Will not merge the multivariate HMMs because their IDs differ.")
        }
    }
        
    ## Check if posteriors are present everywhere
    post.present <- Reduce('|', unlist(lapply(multi.hmm.list, function(x) { !is.null(x$bins$posteriors) })))

    ## Variables
    num.models <- length(multi.hmm.list)
    
    ## Construct list of bins and segments
    message("Concatenating HMMs ...", appendLF=FALSE); ptm <- proc.time()
    bins <- list()    # do not use GRangesList() because it copies the whole list each time an element is added
    segments <- list()
    peaks <- list()
    for (i1 in 1:num.models) {
        hmm <- multi.hmm.list[[1]]    # select always first because we remove it at the end of the loop
        if (!post.present) {
            hmm$bins$posteriors <- NULL
            hmm$segments$mean.posteriors <- NULL
        }
        bins[[i1]] <- hmm$bins
        segments[[i1]] <- hmm$segments
        peaks[[i1]] <- hmm$peaks
        # Remove current HMM to save memory
        if (i1 < num.models) remove(hmm)    # remove it because otherwise R will make a copy when we NULL the underlying reference (multi.hmm.list[[1]])
        multi.hmm.list[[1]] <- NULL
    }
    stopTimedMessage(ptm)

    ## Merge the list
    ptm <- startTimedMessage("Merging ...")
    bins <- do.call('c', bins)    # this can be too memory intensive if posteriors are present
    segments <- do.call('c', segments)
    if (length(peaks) == 1) { # only one chromosome
        peaks.merged <- peaks[[1]]
    } else if (length(peaks) > 1) {
        peaks.merged <- list()
        for (i1 in 1:length(peaks[[1]])) {
            peaks.merged[[names(peaks[[1]])[i1]]] <- do.call('c', lapply(peaks, '[[', i1))
        }
    }
    stopTimedMessage(ptm)

    ## Reassign
    multi.hmm <- hmm
    multi.hmm$bins <- bins
    multi.hmm$segments <- segments
    multi.hmm$peaks <- peaks.merged

    ## Weights
    ptm <- startTimedMessage("Calculating weights ...")
    multi.hmm$weights <- table(multi.hmm$bins$state) / length(multi.hmm$bins)
    stopTimedMessage(ptm)

    if (is.null(filename)) {
        return(multi.hmm)
    } else {
        ptm <- startTimedMessage("Writing to file ",filename," ...")
        save(multi.hmm, file=filename)
        stopTimedMessage(ptm)
    }

    return(NULL)
}
