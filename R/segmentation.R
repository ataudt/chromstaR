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
    df <- as.data.frame(bins)
    ind.readcols <- grep('^counts', names(df))
    ind.postcols <- grep('^posteriors', names(df))
    ind.peakcols <- grep('^peakScores', names(df))
    ind.widthcol <- grep('width', names(df))
    ind.scorecol <- grep('differential.score', names(df))
    red.df <- suppressMessages(collapseBins(df, column2collapseBy=column2collapseBy, columns2average=c(ind.scorecol), columns2drop=c(ind.readcols, ind.widthcol, ind.postcols)))
    names(red.df) <- sub('^mean.','', names(red.df))
    red.peakscores <- as.matrix(red.df[grep('^peakScores', names(red.df))])
    red.df[grep('^peakScores', names(red.df))] <- NULL
    segments <- methods::as(red.df, 'GRanges')
    segments$peakScores <- red.peakscores
    colnames(segments$peakScores) <- colnames(bins$peakScores)
    ## Reorder properly to match the order in bins
    diffscore.temp <- segments$differential.score
    segments$differential.score <- NULL
    segments$differential.score <- diffscore.temp
    segments <- keepSeqlevels(segments, seqlevels(bins))
    seqlengths(segments) <- seqlengths(bins)[seqlevels(segments)]
    stopTimedMessage(ptm)

    return(segments)

}
