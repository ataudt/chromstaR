#' Multivariate segmentation
#'
#' Make segmentation from bins for a \code{\link{multiHMM}} object.
#'
#' @param bins A \code{\link[GenomicRanges]{GRanges}} with binned read counts.
#' @inheritParams collapseBins
#' @return A \code{\link[GenomicRanges]{GRanges}} with segmented regions.
#'
multivariateSegmentation <- function(bins, column2collapseBy='state') {

    ptm <- startTimedMessage("Making segmentation ...")
    ## Intercept case when only one column in posteriors
    if (!is.null(bins$posteriors)) {
        if (ncol(bins$posteriors) == 1) {
            colnames(bins$posteriors) <- paste0('posteriors.', colnames(bins$posteriors))
        }
    }
    df <- as.data.frame(bins)
    if (dim(bins$counts.rpkm)[2] == 1) { # only one experiment
        names2change <- setdiff(names(df)[6:ncol(df)], c('state', 'differential.score', 'combination'))
        names(df)[names(df) %in% names2change] <- gsub('-', '.', paste(setdiff(names(mcols(bins)), c('state', 'differential.score', 'combination')), dimnames(bins$counts.rpkm)[[2]], sep='.'))
      
    }
    ind.readcols <- grep('^counts', names(df))
    ind.postcols <- grep('^posteriors', names(df))
    ind.maxpostcols <- grep('^maxPostInPeak', names(df))
    ind.scorecol <- grep('^differential.score', names(df))
    ind.widthcol <- grep('width', names(df))
    red.df <- suppressMessages(collapseBins(df, column2collapseBy=column2collapseBy, columns2getMax=c(ind.postcols), columns2drop=c(ind.readcols, ind.widthcol, ind.postcols, ind.scorecol)))
    names(red.df) <- sub('^mean\\.','', names(red.df))
    names(red.df) <- sub('^max\\.','', names(red.df))
    red.maxPostInPeak <- NULL
    if (!is.null(bins$posteriors)) {
        red.maxPostInPeak <- as.matrix(red.df[grep('^posteriors', names(red.df))])
        colnames(red.maxPostInPeak) <- colnames(bins$posteriors)
    } else if (!is.null(bins$maxPostInPeak)) {
        red.maxPostInPeak <- as.matrix(red.df[grep('^maxPostInPeak', names(red.df))])
        colnames(red.maxPostInPeak) <- colnames(bins$maxPostInPeak)
    }
    red.df[grep('^posteriors', names(red.df))] <- NULL
    red.df[grep('^maxPostInPeak', names(red.df))] <- NULL
    segments <- methods::as(red.df, 'GRanges')
    ## Null the maxPostInPeak entries outside of peaks
    if (!is.null(red.maxPostInPeak)) {
        segments$maxPostInPeak <- red.maxPostInPeak * dec2bin(segments$state, colnames = colnames(bins$posteriors))
    }
    
    ## Reorder the seqlevels to match the order in bins
    segments <- keepSeqlevels(segments, seqlevels(bins))
    seqlengths(segments) <- seqlengths(bins)[seqlevels(segments)]
    stopTimedMessage(ptm)

    return(segments)

}
