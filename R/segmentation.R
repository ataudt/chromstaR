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
    ind.widthcol <- grep('width', names(df))
    ind.scorecol <- grep('differential.score', names(df))
    red.df <- suppressMessages(collapseBins(df, column2collapseBy=column2collapseBy, columns2average=c(ind.scorecol), columns2drop=c(ind.readcols, ind.widthcol, ind.postcols)))
    names(red.df) <- sub('^mean.','', names(red.df))
    segments <- as(red.df, 'GRanges')
    segments <- keepSeqlevels(segments, seqlevels(bins))
    seqlengths(segments) <- seqlengths(bins)
    stopTimedMessage(ptm)

    return(segments)

}
