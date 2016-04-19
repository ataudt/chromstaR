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
#' @param averages Whether or not averaged posteriors should appear in the segmentation.
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
#'hmm <- callPeaksUnivariate(binned, ID='example_H3K27me3', max.time=60, eps=1,
#'                           keep.posteriors=TRUE)
#'## Compare fits with different post.cutoffs
#'plot(changePostCutoff(hmm, post.cutoff=0.01), type='histogram') + ylim(0,0.25)
#'plot(hmm, type='histogram') + ylim(0,0.25)
#'plot(changePostCutoff(hmm, post.cutoff=0.99), type='histogram') + ylim(0,0.25)
#'
changePostCutoff <- function(model, post.cutoff=0.5, separate.zeroinflation=TRUE, averages=TRUE) {

    post.cutoff <- suppressWarnings( as.numeric(post.cutoff) )
    if (!is.numeric(post.cutoff) | is.na(post.cutoff)) {
        warning("Not changing post.cutoff because given post.cutoff is not numeric.")
        return(model)
    } else if (post.cutoff < 0 | post.cutoff > 1) {
        warning("post.cutoff has to be inside the interval [0,1]. Nothing done.")
        return(model)
    }

    ## post.cutoff threshold
    threshold <- post.cutoff
    model$post.cutoff <- post.cutoff

    ### Univariate HMM ###
    if (is(model,class.univariate.hmm)) {
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
        if (averages==TRUE) {
            red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('score'), columns2drop=c('width',grep('posteriors', names(df), value=TRUE), 'counts')))
            red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'], score=red.df[,'mean.score'])
        } else {
            red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width',grep('posteriors', names(df), value=TRUE),'counts','score')))
            red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'])
        }    
        model$segments <- red.gr
        seqlengths(model$segments) <- seqlengths(model$bins)
        stopTimedMessage(ptm)
#         ## Redo weights
#         model$weights <- table(model$bins$state) / length(model$bins)

    ### Multivariate HMM ###
    } else if (is(model,class.multivariate.hmm)) {
        if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksMultivariate' again with option 'keep.posteriors' set to TRUE.")
        ## Calculate states
        ptm <- startTimedMessage("Calculating states from posteriors ...")
        states <- factor(bin2dec(model$bins$posteriors >= threshold), levels=levels(model$bins$state))
        model$bins$state <- states
        ## Combinations
        if (!is.null(model$bins$combination)) {
            mapping <- model$mapping
#             mapping <- levels(model$bins$combination)
#             names(mapping) <- levels(model$bins$state)
            model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
        }
        stopTimedMessage(ptm)
        ## Redo segmentation
        ptm <- startTimedMessage("Making segmentation ...")
        df <- as.data.frame(model$bins)
        ind.readcols <- grep('^counts', names(df))
        ind.postcols <- grep('^posteriors', names(df))
        ind.widthcol <- grep('width', names(df))
        ind.scorecol <- grep('score', names(df))
        if (averages==TRUE) {
            red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c(ind.postcols,ind.scorecol), columns2drop=c(ind.readcols, ind.widthcol)))
            mean.posteriors <- matrix(unlist(red.df[,grepl('^mean.posteriors',names(red.df))]), ncol=length(model$IDs))
            colnames(mean.posteriors) <- model$IDs
            mean.score <- red.df[,grepl('^mean.score', names(red.df))]
            red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'], score=mean.score)
            red.gr$mean.posteriors <- mean.posteriors
            if (!is.null(model$bins$combination)) {
                red.gr$combination <- red.df[,'combination']
            }
        } else {
            red.df <- suppressMessages(collapseBins(df[,-c(4, ind.readcols, ind.postcols)], column2collapseBy='state'))
            red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'])
            if (!is.null(model$bins$combination)) {
                red.gr$combination <- red.df[,'combination']
            }
        }
        model$segments <- red.gr
        seqlengths(model$segments) <- seqlengths(model$bins)
        stopTimedMessage(ptm)
#         ## Redo weights
#         model$weights <- table(model$bins$state) / length(model$bins)
    } else {
        stop("Supply either a univariate or multivariate chromstaR model")
    }
    ## Return model
    return(model)
    
}

