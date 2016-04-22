#' Calculate differential score
#'
#' Calculate differential score from posteriors.
#'
#' @param posteriors A matrix with posterior probabilities.
#' @param info An \code{\link{experiment.table}} with additional column 'ID'.
#' @return A vector with differential score.
#'
differentialScore <- function(posteriors, info) {

    if (is.null(posteriors) | is.null(info)) {
        return(NULL)
    }
    info.rep <- split(info, as.character(info$replicate))[[1]]
    info.mark <- split(info.rep, as.character(info.rep$mark))
    differential.scores <- list()
    for (mark in names(info.mark)) {
        if (nrow(info.mark[[mark]]) > 0) {
            min.post <- do.call(pmin, as.list(as.data.frame(posteriors[,info.mark[[mark]]$ID])))
            max.post <- do.call(pmax, as.list(as.data.frame(posteriors[,info.mark[[mark]]$ID])))
            differential.scores[[mark]] <- max.post - min.post
        }
    }
    differential.score <- do.call(pmax, differential.scores)
    return(differential.score)
}

