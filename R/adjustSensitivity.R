#' Adjust sensitivity of peak detection
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}}, \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object with a cutoff on the maximum-posterior within each peak. Lower values of \code{sensitivity} mean less sensitive and more precise peak calls. Remaining peaks are kept intact, as opposed to function \code{\link{changePostCutoff}}, where broad peaks are fragmented.
#'
#' Each peak has a maximum-posterior (maxPostInPeak, between 0 and 1) associated. The sensitivity is adjusted with a simple cutoff on 1-maxPostInPeak, e.g. for \code{sensitivity = 0.01} only peaks with a maximum-posterior \code{>= 0.99} will be selected.
#' 
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param sensitivity A vector of values between 0 and 1 for each column in \code{model$bins$posteriors}. If only one value is given, it will be reused for all columns. Values close to 0 will yield more stringent peak calls with lower false positive but higher false negative rate (i.e. more precise but less sensitive).
#' @param invert Select peaks below (\code{FALSE}) or above (\code{TRUE}) the given \code{sensitivity}. This is useful to select low confidence peaks.
#' @return The input object is returned with adjusted peak calls.
#' @export
#' @seealso \code{\link{changePostCutoff}}
#' @examples
#'## Get an example BAM file
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
#'                   stepsizes=500, chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
#'## Compare fits with different fdrs
#'plotHistogram(hmm) + ylim(0,0.25) + ylim(0,0.3)
#'plotHistogram(adjustSensitivity(hmm, sensitivity=0.01)) + ylim(0,0.3)
#'plotHistogram(adjustSensitivity(hmm, sensitivity=1e-12)) + ylim(0,0.3)
#'
adjustSensitivity <- function(model, sensitivity=0.01, invert=FALSE) {

    if (!is.numeric(sensitivity)) {
        warning("'sensitivity' is not numeric. Nothing done.")
        return(model)
    }
    if (is(model, class.univariate.hmm)) {
        model.new <- adjustSensitivity.univariate(model=model, sensitivity=sensitivity, invert=invert)
    } else if (is(model, class.multivariate.hmm) | is(model, class.combined.multivariate.hmm)) {
        model.new <- adjustSensitivity.multivariate(model=model, sensitivity=sensitivity, invert=invert)
    }
    return(model.new)
    
}


adjustSensitivity.multivariate <- function(model, sensitivity, invert=FALSE) {
  
    ## Make sensitivity vector
    if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing.")
    numcol <- ncol(model$bins$posteriors)
    if (length(sensitivity) == 1) {
        sensitivity <- rep(sensitivity, numcol)
    }
    if (length(sensitivity) != numcol) {
        stop("Need ", numcol, " values in 'sensitivity' but only ", length(sensitivity), " are provided.")
    }
    
    for (sensitivity.i in sensitivity) {
        if (sensitivity.i < 0 | sensitivity.i > 1) {
            stop("Values for 'sensitivity' need to be inside the interval [0,1].")
        }
    }
    names(sensitivity) <- colnames(model$bins$posteriors)

    # Check if replicates have the same sensitivity value
    reps <- sub('-rep.*', '', model$info$ID)
    if (any(sapply(split(sensitivity, reps), function(x) { Reduce('&', x==x[1]) }) == FALSE)) {
        stop("Replicates must have the same sensitivity value.")
    }
    
    ## sensitivity threshold
    threshold <- 1-sensitivity

    ### Multivariate HMM ###
    if (is(model, class.multivariate.hmm)) {
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from maximum-posterior in each peak ...")
        p <- getMaxPostInPeaks(model$bins$state, model$bins$posteriors)
        p.thresholded <- matrix(FALSE, ncol=ncol(p), nrow=nrow(p))
        for (icol in 1:ncol(p)) {
            if (!invert) {
                p.thresholded[,icol] <- p[,icol] >= threshold[icol]
            } else {
                p.thresholded[,icol] <- p[,icol] < threshold[icol] & p[,icol] > 0
            }
        }
        states <- factor(bin2dec(p.thresholded), levels=levels(model$bins$state))
        model$bins$state <- states
        ## Combinations
        if (!is.null(model$bins$combination)) {
            mapping <- model$mapping
            model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
        }
        stopTimedMessage(ptm)
        ## Redo peakScores
        ptm <- startTimedMessage("Re-estimating peakScores ...")
        p <- getMaxPostInPeaks(model$bins$state, model$bins$posteriors)
        model$bins$peakScores <- calculatePeakScores(p)
        rm(p)
        stopTimedMessage(ptm)
        ## Redo segmentation
        model$segments <- multivariateSegmentation(model$bins, column2collapseBy='state')
        ## Redo peaks
        ptm <- startTimedMessage("Recalculating peaks ...")
        model$peaks <- list()
        for (i1 in 1:ncol(model$segments$peakScores)) {
            mask <- model$segments$peakScores[,i1] > 0
            peaks <- model$segments[mask]
            mcols(peaks) <- NULL
            peaks$peakScores <- model$segments$peakScores[mask,i1]
            model$peaks[[i1]] <- peaks
        }
        names(model$peaks) <- colnames(model$segments$peakScores)
        stopTimedMessage(ptm)
        
    } else if (is(model, class.combined.multivariate.hmm)) {
        ptm <- startTimedMessage("Calculating states from maximum-posterior in each peak ...")
        mapping.df <- stateBrewer(model$info[,setdiff(names(model$info), 'ID')], mode='full')
        mapping <- mapping.df$combination
        names(mapping) <- mapping.df$state
        p <- getMaxPostInPeaks(model$bins$state, model$bins$posteriors)
        p.thresholded <- matrix(FALSE, ncol=ncol(p), nrow=nrow(p))
        for (icol in 1:ncol(p)) {
            if (!invert) {
                p.thresholded[,icol] <- p[,icol] >= threshold[icol]
            } else {
                p.thresholded[,icol] <- p[,icol] < threshold[icol] & p[,icol] > 0
            }
        }
        states <- factor(bin2dec(p.thresholded), levels=mapping.df$state)
        stopTimedMessage(ptm)
        ## Make fake multiHMM to do the segmentation with "combineMultivariates"
        multiHMM <- list()
        class(multiHMM) <- class.multivariate.hmm
        multiHMM$info <- model$info
        multiHMM$bins <- model$bins
        multiHMM$bins$combination <- mapping[as.character(states)]
        multiHMM$bins$state <- factor(states)
        multiHMM$bins$posteriors <- model$bins$posteriors
        ## Redo peakScores
        ptm <- startTimedMessage("Re-estimating peakScores ...")
        p <- getMaxPostInPeaks(multiHMM$bins$state, multiHMM$bins$posteriors)
        multiHMM$bins$peakScores <- calculatePeakScores(p)
        rm(p)
        stopTimedMessage(ptm)
        multiHMM$mapping <- mapping
        multiHMM$peaks <- model$peaks
        model <- combineMultivariates(list(multiHMM), mode='full')
        ## Redo peaks
        ptm <- startTimedMessage("Recalculating peaks ...")
        model$peaks <- list()
        for (i1 in 1:ncol(model$segments$peakScores)) {
            mask <- model$segments$peakScores[,i1] > 0
            peaks <- model$segments[mask]
            mcols(peaks) <- NULL
            peaks$peakScores <- model$segments$peakScores[mask,i1]
            model$peaks[[i1]] <- peaks
        }
        names(model$peaks) <- colnames(model$segments$peakScores)
        stopTimedMessage(ptm)
    } else {
        stop("Supply either a uniHMM, multiHMM or combinedMultiHMM object.")
    }
    ## Return model
    model$sensitivity <- sensitivity
    return(model)
    

    
}


adjustSensitivity.univariate <- function(model, sensitivity, invert=FALSE) {
  
    sensitivity <- suppressWarnings( as.numeric(sensitivity) )
    if (!is.numeric(sensitivity) | is.na(sensitivity)) {
        warning("Not changing false discovery rate because given 'sensitivity' is not numeric.")
        return(model)
    } else if (sensitivity < 0 | sensitivity > 1) {
        warning("sensitivity has to be inside the interval [0,1]. Nothing done.")
        return(model)
    }

    ## sensitivity threshold
    threshold <- 1-sensitivity

    if (is.null(model$bins$posterior.modified)) stop("Cannot recalculate states because column 'posterior.modified' is missing.")
    ## Calculate states
    ptm <- startTimedMessage("Calculating states from maximum-posterior in each peak ...")
    model$bins$maxPostInPeak <- getMaxPostInPeaks.univariate(model$bins$state, model$bins$posterior.modified)
    states <- model$bins$state
    states[model$bins$state == 'modified'] <- 'unmodified'
    if (!invert) {
        states[ model$bins$maxPostInPeak >= threshold ] <- 'modified'
    } else {
        states[ model$bins$maxPostInPeak < threshold & model$bins$maxPostInPeak > 0 ] <- 'modified'
    }
    model$bins$state <- states
    model$bins$maxPostInPeak <- NULL
    stopTimedMessage(ptm)
    ## Redo peakScores
    ptm <- startTimedMessage("Re-estimating peakScores ...")
    model$bins$maxPostInPeak <- getMaxPostInPeaks.univariate(model$bins$state, model$bins$posterior.modified)
    model$bins$peakScores <- calculatePeakScores.univariate(model$bins$maxPostInPeak)
    model$bins$maxPostInPeak <- NULL
    stopTimedMessage(ptm)
    ## Redo segmentation
    ptm <- startTimedMessage("Making segmentation ...")
    gr <- model$bins
    df <- as.data.frame(gr)
    red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('score'), columns2drop=c('width',grep('posteriors', names(df), value=TRUE), 'counts', 'counts.rpkm', 'posterior.modified')))
    model$segments <- methods::as(red.df, 'GRanges')
    seqlengths(model$segments) <- seqlengths(model$bins)[seqlevels(model$segments)]
    ## Redo peaks
    model$peaks <- model$segments[model$segments$state == 'modified']
    model$peaks$state <- NULL
    model$segments <- NULL
    ## Redo weights
    model$weights <- table(model$bins$state) / length(model$bins)
    stopTimedMessage(ptm)

    ## Return model
    model$sensitivity <- sensitivity
    return(model)
  
}


#' Change the posterior cutoff of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}}, \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object with the given posterior cutoff.
#'
#' Posterior probabilities are between 0 and 1. Peaks are called if the posteriors for a state (univariate) or sample (multivariate) are >= \code{post.cutoff}.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param post.cutoff A vector of posterior cutoff values between 0 and 1 the same length as \code{ncol(model$bins$posteriors)}. If only one value is given, it will be reused for all columns. Values close to 1 will yield more stringent peak calls with lower false positive but higher false negative rate.
#' @return The input object is returned with adjusted peak calls.
#' @export
#' @seealso \code{\link{adjustSensitivity}}
#' @examples
#'## Get an example BAM file
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
#'                   stepsizes=500, chromosomes='chr12')
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
        ## Redo peakScores
        ptm <- startTimedMessage("Re-estimating peakScores ...")
        p <- getMaxPostInPeaks(model$bins$state, model$bins$posteriors)
        model$bins$peakScores <- calculatePeakScores(p)
        rm(p)
        stopTimedMessage(ptm)
        ## Redo segmentation
        model$segments <- multivariateSegmentation(model$bins, column2collapseBy='state')
        ## Redo peaks
        ptm <- startTimedMessage("Recalculating peaks ...")
        model$peaks <- list()
        for (i1 in 1:ncol(model$segments$peakScores)) {
            mask <- model$segments$peakScores[,i1] > 0
            peaks <- model$segments[mask]
            mcols(peaks) <- NULL
            peaks$peakScores <- model$segments$peakScores[mask,i1]
            model$peaks[[i1]] <- peaks
        }
        names(model$peaks) <- colnames(model$segments$peakScores)
        stopTimedMessage(ptm)
    } else if (is(model, class.combined.multivariate.hmm)) {
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        mapping.df <- stateBrewer(model$info[,setdiff(names(model$info), 'ID')], mode='full')
        mapping <- mapping.df$combination
        names(mapping) <- mapping.df$state
        post <- model$bins$posteriors
        post.thresholded <- matrix(FALSE, ncol=ncol(post), nrow=nrow(post))
        for (icol in 1:ncol(post)) {
            post.thresholded[,icol] <- post[,icol] >= threshold[icol]
        }
        states <- factor(bin2dec(post.thresholded), levels=mapping.df$state)
        stopTimedMessage(ptm)
        ## Make fake multiHMM
        multiHMM <- list()
        class(multiHMM) <- class.multivariate.hmm
        multiHMM$info <- model$info
        multiHMM$bins <- model$bins
        multiHMM$bins$combination <- mapping[as.character(states)]
        multiHMM$bins$state <- factor(states)
        multiHMM$bins$posteriors <- post
        ## Redo peakScores
        ptm <- startTimedMessage("Re-estimating peakScores ...")
        p <- getMaxPostInPeaks(multiHMM$bins$state, multiHMM$bins$posteriors)
        multiHMM$bins$peakScores <- calculatePeakScores(p)
        rm(p)
        stopTimedMessage(ptm)
        multiHMM$mapping <- mapping
        multiHMM$peaks <- model$peaks
        model <- combineMultivariates(list(multiHMM), mode='full')
        ## Redo peaks
        ptm <- startTimedMessage("Recalculating peaks ...")
        model$peaks <- list()
        for (i1 in 1:ncol(model$segments$peakScores)) {
            mask <- model$segments$peakScores[,i1] > 0
            peaks <- model$segments[mask]
            mcols(peaks) <- NULL
            peaks$peakScores <- model$segments$peakScores[mask,i1]
            model$peaks[[i1]] <- peaks
        }
        names(model$peaks) <- colnames(model$segments$peakScores)
        stopTimedMessage(ptm)
    } else {
        stop("Supply either a uniHMM, multiHMM or combinedMultiHMM object.")
    }
    
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
    ## Redo peakScores
    ptm <- startTimedMessage("Re-estimating peakScores ...")
    model$bins$maxPostInPeak <- getMaxPostInPeaks.univariate(model$bins$state, model$bins$posterior.modified)
    model$bins$peakScores <- calculatePeakScores.univariate(model$bins$maxPostInPeak)
    model$bins$maxPostInPeak <- NULL
    stopTimedMessage(ptm)
    ## Redo segmentation
    ptm <- startTimedMessage("Making segmentation ...")
    gr <- model$bins
    df <- as.data.frame(gr)
    red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width', 'posterior.modified', grep('posteriors', names(df), value=TRUE), 'counts', 'counts.rpkm', 'posterior.modified')))
    segments <- methods::as(red.df, 'GRanges')
    model$peaks <- segments[segments$state == 'modified']
    model$peaks$state <- NULL
    seqlengths(model$peaks) <- seqlengths(model$bins)[seqlevels(model$peaks)]
    stopTimedMessage(ptm)

    ## Return model
    return(model)
  
}
