#' Remove condition from model
#' 
#' Remove a condition from a \code{\link{combinedMultiHMM}} object.
#' 
#' @param model A \code{\link{combinedMultiHMM}} object or file which contains such an object.
#' @param conditions A character vector with the condition(s) to be removed.
#' @return The input \code{\link{combinedMultiHMM}} object with specified conditions removed.
#' @export
#' @examples 
#'## Get an example HMM
#'file <- system.file("data","combined_mode-differential.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'
#'## Print available conditions
#'print(unique(model$info$condition))
#'
#'## Remove condition SHR
#'new.model <- removeCondition(model, conditions='SHR')
#'
removeCondition <- function(model, conditions) {
  
    model <- loadHmmsFromFiles(model, check.class=class.combined.multivariate.hmm)[[1]]
    
    ### Check input ###
    cond.setdiff <- setdiff(conditions, model$info$condition)
    if (length(cond.setdiff) == length(conditions)) {
        warning("Conditions ", paste0(cond.setdiff, collapse=', '), " could not be found. Doing nothing.")
        return(model)
    }
    if (length(cond.setdiff) > 0) {
        warning("Conditions ", paste0(cond.setdiff, collapse=', '), " could not be found.")
    }
    conditions <- setdiff(conditions, cond.setdiff)
    
    ### Info ###
    info <- model$info
    mask <- !info$condition %in% conditions
    info <- info[mask,]
    model$info <- info
    
    ### Bins ###
    # State
    if (!is.null(model$bins$state)) {
        ptm <- startTimedMessage("Recalculating states ...")
        states <- model$bins$state
        binstates <- dec2bin(states, colnames=model$info$ID)
        removeconds <- paste0(paste0('-', conditions, '-'), collapse='|')
        binstates <- binstates[, grep(removeconds, colnames(binstates), invert=TRUE)]
        state <- factor(bin2dec(binstates))
        model$bins$state <- state
        stopTimedMessage(ptm)
    }
    # Combinations
    ptm <- startTimedMessage("Selecting conditions ...")
    removeconds <- paste0('combination.', conditions)
    mcols(model$bins)[removeconds] <- NULL
    # Counts
    if (!is.null(model$bins$counts)) {
        counts <- model$bins$counts
        removeconds <- paste0(paste0('-', conditions, '-'), collapse='|')
        keepconds <- grep(removeconds, colnames(counts), invert=TRUE, value=TRUE)
        counts <- counts[,keepconds]
        if (class(counts) != "matrix") {
            counts <- matrix(counts, ncol=1, dimnames=list(NULL, keepconds))
        }
        model$bins$counts <- counts
    }
    if (!is.null(model$bins$posteriors)) {
        posteriors <- model$bins$posteriors
        removeconds <- paste0(paste0('-', conditions, '-'), collapse='|')
        keepconds <- grep(removeconds, colnames(posteriors), invert=TRUE, value=TRUE)
        posteriors <- posteriors[,keepconds]
        if (class(posteriors) != "matrix") {
            posteriors <- matrix(posteriors, ncol=1, dimnames=list(NULL, keepconds))
        }
        model$bins$posteriors <- posteriors
    }
    if (!is.null(model$bins$peakScores)) {
        peakScores <- model$bins$peakScores
        removeconds <- paste0(paste0('-', conditions, '-'), collapse='|')
        keepconds <- grep(removeconds, colnames(peakScores), invert=TRUE, value=TRUE)
        peakScores <- peakScores[,keepconds]
        if (class(peakScores) != "matrix") {
            peakScores <- matrix(peakScores, ncol=1, dimnames=list(NULL, keepconds))
        }
        model$bins$peakScores <- peakScores
    }
    stopTimedMessage(ptm)
    # Redo differential score
    if (!is.null(model$bins$differential.score)) {
        ptm <- startTimedMessage("Differential score ...")
        model$bins$differential.score <- differentialScoreSum(model$bins$peakScores, model$info)
        stopTimedMessage(ptm)
    }
    # Redo transition frequencies and groups
    ptm <- startTimedMessage("Transition frequencies ...")
    trans <- suppressMessages( transitionFrequencies(combined.hmm = model) )
    model$bins$transition.group <- trans$per.bin$group
    model$frequencies <- trans$table
    stopTimedMessage(ptm)
    
    bins <- model$bins
    combs <- getCombinations(model$bins)
    names(combs) <- sub('combination.', '', names(combs))
    ### Redo the segmentation for all conditions combined ###
    ptm <- startTimedMessage("Redoing segmentation for all conditions combined ...")
    segments <- suppressMessages( multivariateSegmentation(bins, column2collapseBy='state') )
    names(mcols(segments)) <- setdiff(names(mcols(bins)), c('posteriors','counts'))
    model$segments <- segments
    stopTimedMessage(ptm)
    
    ### Redo the segmentation for each condition separately ###
    ptm <- startTimedMessage("Redoing segmentation for each condition separately ...")
    segments.per.condition <- list()
    for (cond in names(combs)) {
        bins.cond <- bins
        mcols(bins.cond) <- mcols(bins)[paste0('combination.',cond)]
        df <- as.data.frame(bins.cond)
        names(df)[6] <- cond
        segments.cond <- suppressMessages( collapseBins(df, column2collapseBy=cond, columns2drop=c('width', grep('posteriors', names(df), value=TRUE))) )
        segments.cond <- methods::as(segments.cond, 'GRanges')
        names(mcols(segments.cond)) <- 'combination'
        seqlengths(segments.cond) <- seqlengths(bins)[seqlevels(segments.cond)]
        segments.per.condition[[cond]] <- segments.cond
    }
    model$segments.per.condition <- segments.per.condition
    stopTimedMessage(ptm)
    
    ### Peaks ###
    removeconds <- paste0(paste0('-', conditions, '-'), collapse='|')
    peaks <- model$peaks[grep(removeconds, names(model$peaks), invert=TRUE)]
    model$peaks <- peaks
    
    return(model)
}
