#' chromstaR scores
#'
#' Various scores used in \code{\link{chromstaR}}.
#'
#' @param mat A matrix with posterior probabilities, read counts or any other matrix with these dimensions. Column names must correspond to the ID entries in \code{info}.
#' @param info An \code{\link{experiment.table}} with additional column 'ID'.
#' @param FUN A function to compute the score with.
#' @name scores
#' @return A numeric vector.
#' @author Aaron Taudt
NULL

#' @describeIn scores Maximum differential score. Values are between 0 and 1. A value of 1 means that at least one mark is maximally different between conditions.
differentialScoreMax <- function(mat, info, FUN='-') {

    if (is.null(mat) | is.null(info)) {
        return(NULL)
    }
    differential.scores <- differentialScores(mat, info, FUN)
    differential.score <- do.call(pmax, differential.scores)
    return(differential.score)
}


#' @describeIn scores Additive differential score. Values are between 0 and N, where N is the number of marks. A value around 1 means that approximately 1 mark is different, a value of 2 means that 2 marks are different etc.
differentialScoreSum <- function(mat, info, FUN='-') {

    if (is.null(mat) | is.null(info)) {
        return(NULL)
    }
    differential.scores <- differentialScores(mat, info, FUN)
    differential.score <- Reduce('+', differential.scores)
    return(differential.score)
}


differentialScores <- function(mat, info, FUN='-') {

    FUN <- match.fun(FUN)
    if (is.null(mat) | is.null(info)) {
        return(NULL)
    }
    # info.rep <- split(info, paste0(info$mark, '-', info$condition))
    # for (i1 in 1:length(info.rep)) {
    #     info.rep[[i1]]$replicate <- as.integer(as.factor(as.character(info.rep[[i1]]$replicate)))
    # }
    # info.rep <- do.call(rbind, info.rep)
    # info.rep <- split(info.rep, as.character(info.rep$replicate))[[1]]
    # info.mark <- split(info.rep, as.character(info.rep$mark))
    info.mark <- split(info, as.character(info$mark))
    differential.scores <- list()
    for (mark in names(info.mark)) {
        if (nrow(info.mark[[mark]]) > 1) {
            min.post <- do.call(pmin, as.list(as.data.frame(mat[,info.mark[[mark]]$ID,drop=FALSE])))
            max.post <- do.call(pmax, as.list(as.data.frame(mat[,info.mark[[mark]]$ID,drop=FALSE])))
            differential.scores[[mark]] <- FUN(max.post, min.post)
        } else if (nrow(info.mark[[mark]]) == 1) {
            differential.scores[[mark]] <- rep(0, nrow(mat))
        }
    }
    return(differential.scores)
}
