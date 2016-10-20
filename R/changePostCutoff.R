#' Change the posterior cutoff of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}} or \code{\link{multiHMM}} object with the given posterior cutoff.
#'
#' Posterior probabilities are between 0 and 1. Peaks are called if the posteriors for a state (univariate) or sample (multivariate) are >= \code{post.cutoff}.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param post.cutoff A vector of posterior cutoff values between 0 and 1 the same length as \code{ncol(model$bins$posteriors)}. If only one value is given, it will be reused for all columns. Values close to 1 will yield more stringent peak calls with lower false positive but higher false negative rate.
#' @return The input object is returned with adjusted peak calls.
#' @export
#' @examples
#'## Get an example BAM file
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'# !Keep posteriors to change the post.cutoff later!
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1,
#'                           keep.posteriors=TRUE)
#'## Compare fits with different post.cutoffs
#'plotHistogram(changePostCutoff(hmm, post.cutoff=0.01)) + ylim(0,0.25)
#'plotHistogram(hmm) + ylim(0,0.25)
#'plotHistogram(changePostCutoff(hmm, post.cutoff=0.99)) + ylim(0,0.25)
#'
changePostCutoff <- function(model, post.cutoff=0.5) {

    if (!is.numeric(post.cutoff)) {
        warning("'post.cutoff' is not numeric. Nothing done.")
        return(model)
    }
    if (is(model, class.univariate.hmm)) {
        model.new <- changePostCutoff.univariate(model=model, post.cutoff=post.cutoff)
    } else if (is(model, class.multivariate.hmm) | is(model, class.combined.multivariate.hmm)) {
        model.new <- changePostCutoff.multivariate(model=model, post.cutoff=post.cutoff)
    }
    return(model.new)
    
}


changePostCutoff.multivariate <- function(model, post.cutoff) {
  
    ## Make post.cutoff vector
    if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing.")
    numcol <- ncol(model$bins$posteriors)
    if (length(post.cutoff) == 1) {
        post.cutoff <- rep(post.cutoff, numcol)
    }
    if (length(post.cutoff) != numcol) {
        stop("Need ", numcol, " values in 'post.cutoff' but only ", length(post.cutoff), " are provided.")
    }
    
    for (post.cutoff.i in post.cutoff) {
        if (post.cutoff.i < 0 | post.cutoff.i > 1) {
            stop("Values for 'post.cutoff' need to be inside the interval [0,1].")
        }
    }
    names(post.cutoff) <- colnames(model$bins$posteriors)

    # Check if replicates have the same post.cutoff value
    reps <- sub('-rep.*', '', model$info$ID)
    if (any(sapply(split(post.cutoff, reps), function(x) { Reduce('&', x==x[1]) }) == FALSE)) {
        stop("Replicates must have the same post.cutoff value.")
    }
    
    ## post.cutoff threshold
    threshold <- post.cutoff

    ### Multivariate HMM ###
    if (is(model, class.multivariate.hmm)) {
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        post <- model$bins$posteriors
        post.thresholded <- matrix(FALSE, ncol=ncol(post), nrow=nrow(post))
        for (icol in 1:ncol(post)) {
            post.thresholded[,icol] <- post[,icol] >= threshold[icol]
        }
        states <- factor(bin2dec(post.thresholded), levels=levels(model$bins$state))
        model$bins$state <- states
        ## Combinations
        if (!is.null(model$bins$combination)) {
            mapping <- model$mapping
            model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
        }
        stopTimedMessage(ptm)
        ## Redo segmentation
        model$segments <- multivariateSegmentation(model$bins, column2collapseBy='state')
    } else if (is(model, class.combined.multivariate.hmm)) {
        mapping.df <- stateBrewer(model$info[,setdiff(names(model$info), 'ID')], mode='full')
        mapping <- mapping.df$combination
        names(mapping) <- mapping.df$state
        post <- model$bins$posteriors
        post.thresholded <- matrix(FALSE, ncol=ncol(post), nrow=nrow(post))
        for (icol in 1:ncol(post)) {
            post.thresholded[,icol] <- post[,icol] >= threshold[icol]
        }
        states <- factor(bin2dec(post.thresholded), levels=mapping.df$state)
        ## Make fake multiHMM
        multiHMM <- list()
        class(multiHMM) <- class.multivariate.hmm
        multiHMM$info <- model$info
        multiHMM$bins <- model$bins
        multiHMM$bins$combination <- mapping[as.character(states)]
        multiHMM$bins$state <- factor(states)
        multiHMM$bins$posteriors <- post
        multiHMM$mapping <- mapping
        multiHMM$peaks <- model$peaks
        model <- combineMultivariates(list(multiHMM), mode='full')
    } else {
        stop("Supply either a uniHMM, multiHMM or combinedMultiHMM object.")
    }
    
    ## Redo peaks
    ptm <- startTimedMessage("Recalculating peaks ...")
    binstates <- dec2bin(model$segments$state, colnames=colnames(post))
    segments <- model$segments
    mcols(segments) <- NULL
    for (icol in 1:ncol(binstates)) {
        segments$peakScores <- model$segments$peakScores[,icol]
        segments$state <- binstates[,icol]
        red.df <- suppressMessages( collapseBins(as.data.frame(segments), column2collapseBy = 'state') )
        red.df <- red.df[red.df$state == 1,]
        red.df$state <- NULL
        model$peaks[[icol]] <- methods::as(red.df, 'GRanges')
    }
    stopTimedMessage(ptm)
    
    ## Return model
    model$post.cutoff <- threshold
    return(model)
    

    
}


changePostCutoff.univariate <- function(model, post.cutoff) {
  
    post.cutoff <- suppressWarnings( as.numeric(post.cutoff) )
    if (!is.numeric(post.cutoff) | is.na(post.cutoff)) {
        warning("Not changing posterior cutoff because given 'post.cutoff' is not numeric.")
        return(model)
    } else if (post.cutoff < 0 | post.cutoff > 1) {
        warning("post.cutoff has to be inside the interval [0,1]. Nothing done.")
        return(model)
    }

    ## post.cutoff threshold
    threshold <- post.cutoff
    model$post.cutoff <- post.cutoff

    if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksUnivariate' again with option 'keep.posteriors' set to TRUE.")
    ## Calculate states
    ptm <- startTimedMessage("Calculating states from posteriors ...")
    states <- rep(NA,length(model$bins))
    states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]<=model$bins$posteriors[,1] ] <- 1
    states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]>model$bins$posteriors[,1] ] <- 2
    states[ model$bins$posteriors[,3]>=threshold ] <- 3
    states <- state.labels[states]
    model$bins$state <- states
    stopTimedMessage(ptm)
    ## Redo segmentation
    ptm <- startTimedMessage("Making segmentation ...")
    gr <- model$bins
    df <- as.data.frame(gr)
    red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width', 'posterior.modified', grep('posteriors', names(df), value=TRUE), 'counts')))
    segments <- methods::as(red.df, 'GRanges')
    model$peaks <- segments[segments$state == 'modified']
    model$peaks$state <- NULL
    seqlengths(model$peaks) <- seqlengths(model$bins)[seqlevels(model$peaks)]
    stopTimedMessage(ptm)

    ## Return model
    return(model)
  
}
