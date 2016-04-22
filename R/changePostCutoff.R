#' Change the false discovery rate of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}} or \code{\link{multiHMM}} object with the given posterior cutoff.
#'
#' Posterior probabilities are between 0 and 1. Peaks are called if the posteriors for a state (univariate) or sample (multivariate) are >= \code{post.cutoff}.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param post.cutoff Posterior cutoff. Values close to 1 will yield more stringent peak calls with lower false positive but higher false negative rate.
#' @param separate.zeroinflation Only for \code{\link{uniHMM}} objects: If set to TRUE, state 'zero-inflation' will be treated separately, otherwise it will be merged with state 'unmodified'.
#' @return The input object is returned with adjusted peak calls.
#' @export
#' @examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(bedfile, assembly=rn4_chrominfo, binsize=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'# !Keep posteriors to change the post.cutoff later!
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1,
#'                           keep.posteriors=TRUE)
#'## Compare fits with different post.cutoffs
#'plot(changePostCutoff(hmm, post.cutoff=0.01), type='histogram') + ylim(0,0.25)
#'plot(hmm, type='histogram') + ylim(0,0.25)
#'plot(changePostCutoff(hmm, post.cutoff=0.99), type='histogram') + ylim(0,0.25)
#'
changePostCutoff <- function(model, post.cutoff=0.5, separate.zeroinflation=TRUE) {

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

    ### Univariate HMM ###
    if (is(model, class.univariate.hmm)) {
        if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksUnivariate' again with option 'keep.posteriors' set to TRUE.")
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        states <- rep(NA,length(model$bins))
        if (separate.zeroinflation) {
            states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]<=model$bins$posteriors[,1] ] <- 1
            states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]>model$bins$posteriors[,1] ] <- 2
            states[ model$bins$posteriors[,3]>=threshold ] <- 3
            states <- state.labels[states]
        } else {
            states <- ifelse(model$bins$posteriors[,3]>=threshold, 2, 1)
            states <- state.labels[2:3][states]
        }
        model$bins$state <- states
        stopTimedMessage(ptm)
        ## Redo segmentation
        ptm <- startTimedMessage("Making segmentation ...")
        gr <- model$bins
        df <- as.data.frame(gr)
        red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('score'), columns2drop=c('width',grep('posteriors', names(df), value=TRUE), 'counts')))
        red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'], score=red.df[,'mean.score'])
        model$segments <- red.gr
        seqlengths(model$segments) <- seqlengths(model$bins)
        stopTimedMessage(ptm)
#         ## Redo weights
#         model$weights <- table(model$bins$state) / length(model$bins)

    ### Multivariate HMM ###
    } else if (is(model, class.multivariate.hmm)) {
        if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksMultivariate' again with option 'keep.posteriors' set to TRUE.")
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        states <- factor(bin2dec(model$bins$posteriors >= threshold), levels=levels(model$bins$state))
        model$bins$state <- states
        ## Combinations
        if (!is.null(model$bins$combination)) {
            mapping <- model$mapping
            model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
        }
        stopTimedMessage(ptm)
        ## Redo segmentation
        model$segments <- multivariateSegmentation(model$bins, column2collapseBy='state')
#         ## Redo weights
#         model$weights <- table(model$bins$state) / length(model$bins)
    } else if (is(model, class.combined.multivariate.hmm)) {
        if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksMultivariate' again with option 'keep.posteriors' set to TRUE.")
        mapping.df <- stateBrewer(model$info[,1:5], mode='full')
        mapping <- mapping.df$combination
        names(mapping) <- mapping.df$state
        post <- model$bins$posteriors
        states <- factor(bin2dec(post >= threshold), levels=mapping.df$state)
        ## Make fake multiHMM
        multiHMM <- list()
        class(multiHMM) <- class.multivariate.hmm
        multiHMM$info <- model$info
        multiHMM$bins <- model$bins
        multiHMM$bins$combination <- mapping[as.character(states)]
        multiHMM$bins$state <- factor(states)
        multiHMM$bins$posteriors <- post
        multiHMM$mapping <- mapping
        model <- combineMultivariates(list(multiHMM), mode='full')
        model$post.cutoff <- threshold
    } else {
        stop("Supply either a uniHMM, multiHMM or combinedMultiHMM object.")
    }
    ## Return model
    return(model)
    
}

