#' chromstaR scores
#'
#' Various scores used in \code{\link{chromstaR}}.
#'
#' @param posteriors A matrix with posterior probabilities.
#' @param info An \code{\link{experiment.table}} with additional column 'ID'.
#' @name scores
#' @return A numeric vector.
#' @author Aaron Taudt
NULL

#' @describeIn scores Maximum differential score. Values are between 0 and 1. A value of 1 means that at least one mark is maximally different between conditions.
differentialScoreMax <- function(posteriors, info) {

    if (is.null(posteriors) | is.null(info)) {
        return(NULL)
    }
    differential.scores <- differentialScores(posteriors, info)
    differential.score <- do.call(pmax, differential.scores)
    return(differential.score)
}


#' @describeIn scores Additive differential score. Values are between 0 and N, where N is the number of marks. A value around 1 means that approximately 1 mark is different, a value of 2 means that 2 marks are different etc.
differentialScoreSum <- function(posteriors, info) {

    if (is.null(posteriors) | is.null(info)) {
        return(NULL)
    }
    differential.scores <- differentialScores(posteriors, info)
    differential.score <- rowSums(as.data.frame(differential.scores))
    return(differential.score)
}


differentialScores <- function(posteriors, info) {

    if (is.null(posteriors) | is.null(info)) {
        return(NULL)
    }
    info.rep <- split(info, paste0(info$mark, '-', info$condition))
    for (i1 in 1:length(info.rep)) {
        info.rep[[i1]]$replicate <- as.integer(as.factor(as.character(info.rep[[i1]]$replicate)))
    }
    info.rep <- do.call(rbind, info.rep)
    info.rep <- split(info.rep, as.character(info.rep$replicate))[[1]]
    info.mark <- split(info.rep, as.character(info.rep$mark))
    differential.scores <- list()
    for (mark in names(info.mark)) {
        if (nrow(info.mark[[mark]]) > 1) {
            min.post <- do.call(pmin, as.list(as.data.frame(posteriors[,info.mark[[mark]]$ID])))
            max.post <- do.call(pmax, as.list(as.data.frame(posteriors[,info.mark[[mark]]$ID])))
            differential.scores[[mark]] <- max.post - min.post
        } else if (nrow(info.mark[[mark]]) == 1) {
            differential.scores[[mark]] <- rep(0, nrow(posteriors))
        }
    }
    return(differential.scores)
}

