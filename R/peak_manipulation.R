getPeakScores <- function(bins) {
    
    binstates <- dec2bin(bins$state, colnames = colnames(bins$posteriors))
    rownames(binstates) <- NULL
    peakScores <- matrix(NA, nrow=nrow(binstates), ncol=ncol(binstates), dimnames=list(NULL, colnames(binstates)))
    for (icol in 1:ncol(binstates)) {
        r.bin <- rle(binstates[,icol])
        r <- r.bin
        r$values[r$values == TRUE] <- 1:length(which(r$values==TRUE))
        r$values[r$values == FALSE] <- NA
        peakNumbers <- inverse.rle(r)
        df <- aggregate(bins$posteriors[,icol], by=list(peakNumber=peakNumbers), FUN=max)
        if (class(df$x) == 'list') {
            class(df$x) <- 'numeric'
        }
        r <- r.bin
        r$values[r$values == TRUE] <- df$x
        r$values[r$values == FALSE] <- 0
        peakScores[,icol] <- inverse.rle(r)
    }
    return(peakScores)
    
}


getPeakScore.univariate <- function(bins) {
    
    p <- bins$posterior.modified
    binstates <- bins$state == 'modified'
    r.bin <- rle(binstates)
    r <- r.bin
    r$values[r$values == TRUE] <- 1:length(which(r$values==TRUE))
    r$values[r$values == FALSE] <- NA
    peakNumbers <- inverse.rle(r)
    df <- aggregate(p, by=list(peakNumber=peakNumbers), FUN=max)
    r <- r.bin
    r$values[r$values == TRUE] <- df$x
    r$values[r$values == FALSE] <- 0
    peakScore <- inverse.rle(r)
    return(peakScore)
    
}