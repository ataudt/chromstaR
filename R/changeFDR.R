#' Change the false discovery rate of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}} or \code{\link{multiHMM}} object with the given false discovery rate.
#'
#' Each peak has a peak score associated (between 0 and 1). The false discovery rate is a simple cutoff on 1-peakScore, e.g. for an \code{fdr = 0.01} only peaks with a peak score \code{>= 0.99} (1-fdr) will be selected.
#' 
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param fdr A vector of values between 0 and 1 for each column in \code{model$bins$peakScores}. If only one value is given, it will be reused for all columns. Values close to 0 will yield more stringent peak calls with lower false positive but higher false negative rate.
#' @param invert Select peaks below (\code{FDR}) or above (\code{TRUE}) the given \code{fdr}. This is useful to select low confidence peaks.
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
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
#'## Compare fits with different fdrs
#'plotHistogram(hmm) + ylim(0,0.25) + ylim(0,0.3)
#'plotHistogram(changeFDR(hmm, fdr=0.01)) + ylim(0,0.3)
#'plotHistogram(changeFDR(hmm, fdr=1e-12)) + ylim(0,0.3)
#'
changeFDR <- function(model, fdr=0.01, invert=FALSE) {

    if (!is.numeric(fdr)) {
        warning("'fdr' is not numeric. Nothing done.")
        return(model)
    }
    if (is(model, class.univariate.hmm)) {
        model.new <- changeFDR.univariate(model=model, fdr=fdr, invert=invert)
    } else if (is(model, class.multivariate.hmm) | is(model, class.combined.multivariate.hmm)) {
        model.new <- changeFDR.multivariate(model=model, fdr=fdr, invert=invert)
    }
    return(model.new)
    
}


changeFDR.multivariate <- function(model, fdr, invert=FALSE) {
  
    ## Make fdr vector
    if (is.null(model$bins$peakScores)) stop("Cannot recalculate states because peakScores are missing.")
    numcol <- ncol(model$bins$peakScores)
    if (length(fdr) == 1) {
        fdr <- rep(fdr, numcol)
    }
    if (length(fdr) != numcol) {
        stop("Need ", numcol, " values in 'fdr' but only ", length(fdr), " are provided.")
    }
    
    for (fdr.i in fdr) {
        if (fdr.i < 0 | fdr.i > 1) {
            stop("Values for 'fdr' need to be inside the interval [0,1].")
        }
    }
    names(fdr) <- colnames(model$bins$peakScores)

    # Check if replicates have the same fdr value
    reps <- sub('-rep.*', '', model$info$ID)
    if (any(sapply(split(fdr, reps), function(x) { Reduce('&', x==x[1]) }) == FALSE)) {
        stop("Replicates must have the same fdr value.")
    }
    
    ## fdr threshold
    threshold <- 1-fdr

    ### Multivariate HMM ###
    if (is(model, class.multivariate.hmm)) {
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from peakScores ...")
        p <- model$bins$peakScores
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
        ## Redo segmentation
        model$segments <- multivariateSegmentation(model$bins, column2collapseBy='state')
        ## Redo peaks
        model$peaks <- list()
        for (i1 in 1:ncol(model$segments$peakScores)) {
            mask <- model$segments$peakScores[,i1] > 0
            peaks <- model$segments[mask]
            mcols(peaks) <- NULL
            peaks$peakScores <- model$segments$peakScores[mask,i1]
            model$peaks[[i1]] <- peaks
        }
        names(model$peaks) <- colnames(model$segments$peakScores)
        
    } else if (is(model, class.combined.multivariate.hmm)) {
        mapping.df <- stateBrewer(model$info[,setdiff(names(model$info), 'ID')], mode='full')
        mapping <- mapping.df$combination
        names(mapping) <- mapping.df$state
        p <- model$bins$peakScores
        p.thresholded <- matrix(FALSE, ncol=ncol(p), nrow=nrow(p))
        for (icol in 1:ncol(p)) {
            if (!invert) {
                p.thresholded[,icol] <- p[,icol] >= threshold[icol]
            } else {
                p.thresholded[,icol] <- p[,icol] < threshold[icol] & p[,icol] > 0
            }
        }
        states <- factor(bin2dec(p.thresholded), levels=mapping.df$state)
        ## Make fake multiHMM
        multiHMM <- list()
        class(multiHMM) <- class.multivariate.hmm
        multiHMM$info <- model$info
        multiHMM$bins <- model$bins
        multiHMM$bins$combination <- mapping[as.character(states)]
        multiHMM$bins$state <- factor(states)
        multiHMM$bins$posteriors <- model$bins$posteriors
        multiHMM$bins$peakScores <- p
        multiHMM$mapping <- mapping
        multiHMM$peaks <- model$peaks
        model <- combineMultivariates(list(multiHMM), mode='full')
        ## Apply FDR on peaks ##
        for (i1 in 1:length(model$peaks)) {
            peaks <- model$peaks[[i1]]
            if (!invert) {
                model$peaks[[i1]] <- peaks[peaks$peakScores >= threshold[names(model$peaks)[i1]]]
            } else {
                model$peaks[[i1]] <- peaks[peaks$peakScores < threshold[names(model$peaks)[i1]]]
            }
        }
    } else {
        stop("Supply either a uniHMM, multiHMM or combinedMultiHMM object.")
    }
    ## Return model
    model$fdr <- fdr
    return(model)
    

    
}


changeFDR.univariate <- function(model, fdr, invert=FALSE) {
  
    fdr <- suppressWarnings( as.numeric(fdr) )
    if (!is.numeric(fdr) | is.na(fdr)) {
        warning("Not changing false discovery rate because given 'fdr' is not numeric.")
        return(model)
    } else if (fdr < 0 | fdr > 1) {
        warning("fdr has to be inside the interval [0,1]. Nothing done.")
        return(model)
    }

    ## fdr threshold
    threshold <- 1-fdr

    if (is.null(model$bins$peakScores)) stop("Cannot recalculate states because column 'peakScores' is missing.")
    ## Calculate states
    ptm <- startTimedMessage("Calculating states from peakScores ...")
    states <- model$bins$state
    states[model$bins$state == 'modified'] <- 'unmodified'
    if (!invert) {
        states[ model$bins$peakScores >= threshold ] <- 'modified'
    } else {
        states[ model$bins$peakScores < threshold & model$bins$peakScores > 0 ] <- 'modified'
    }
    model$bins$state <- states
    stopTimedMessage(ptm)
    ## Redo segmentation
    ptm <- startTimedMessage("Making segmentation ...")
    gr <- model$bins
    df <- as.data.frame(gr)
    red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('score'), columns2drop=c('width',grep('posteriors', names(df), value=TRUE), 'counts')))
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
    model$fdr <- fdr
    return(model)
  
}
