#' Add score for differential export
#'
#' Add column 'score' to the 'segments' entry in a \code{\link{multiHMM}}. This score can be used to rank differential states.
#'
#' @author Aaron Taudt
#' @param multi.model A \code{\link{multiHMM}} object or a file that contains such an object.
#' @inheritParams state.brewer
#' @param reorderSegments If \code{TRUE} segments will be reordered by score.
#' @return The input \code{\link{multiHMM}} with additional column 'score' in \code{multi.model$segments}
addScore <- function(multi.model, tracks2compare, reorderSegments=TRUE) {
    
    multi.model <- loadHmmsFromFiles(multi.model, check.class=class.multivariate.hmm)[[1]]
    segments <- multi.model$segments
    binstates.all <- dec2bin(segments$state, colnames=multi.model$IDs)
    unique.tracks <- unique(tracks2compare)
    ## Go through the different tracks2compare
    score <- list()
    for (track in unique.tracks) {
        binstates <- binstates.all[,track %in% tracks2compare]
        nbinstates <- apply(!binstates, 2, as.integer)
        binstates <- apply(binstates, 2, as.integer)
        ## Variance of peaks to non-peaks of maximum posteriors per segment
        score[[track]] <- apply(cbind(apply(binstates*segments$mean.posteriors, 1, max), apply(nbinstates*segments$mean.posteriors, 1, max)), 1, var)
    }
    score <- do.call(cbind, score)
    ## Add up score between tracks2compare
    segments$score <- apply(score, 1, sum)
    ## Reorder segments
    if (reorderSegments) {
        segments <- segments[rev(order(segments$score))]
    }

    multi.model$segments <- segments
    return(multi.model)
    
}
