getMaxPostInPeaks <- function(states, posteriors) {
    
    binstates <- dec2bin(states, colnames = colnames(posteriors))
    rownames(binstates) <- NULL
    peakScores <- matrix(NA, nrow=nrow(binstates), ncol=ncol(binstates), dimnames=list(NULL, colnames(binstates)))
    for (icol in 1:ncol(binstates)) {
        r.bin <- rle(binstates[,icol])
        r <- r.bin
        r$values[r$values == TRUE] <- 1:length(which(r$values==TRUE))
        r$values[r$values == FALSE] <- NA
        peakNumbers <- inverse.rle(r)
        df <- aggregate(posteriors[,icol], by=list(peakNumber=peakNumbers), FUN=max)
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


getMaxPostInPeaks.univariate <- function(states, posteriors) {
    
    binstates <- states == 'modified'
    r.bin <- rle(binstates)
    r <- r.bin
    r$values[r$values == TRUE] <- 1:length(which(r$values==TRUE))
    r$values[r$values == FALSE] <- NA
    peakNumbers <- inverse.rle(r)
    df <- aggregate(posteriors, by=list(peakNumber=peakNumbers), FUN=max)
    r <- r.bin
    r$values[r$values == TRUE] <- df$x
    r$values[r$values == FALSE] <- 0
    peakScore <- inverse.rle(r)
    return(peakScore)
    
}


calculatePeakScores <- function(maxPostInPeak) {
    peakScores <- maxPostInPeak
    for (i1 in 1:ncol(maxPostInPeak)) {
        mask <- maxPostInPeak[,i1] > 0
        peakScores[mask,i1] <- stats::ecdf(maxPostInPeak[mask,i1])(maxPostInPeak[mask,i1])*1000
    }
    return(peakScores)
}


calculatePeakScores.univariate <- function(maxPostInPeak) {
    peakScores <- maxPostInPeak
    mask <- maxPostInPeak > 0
    peakScores[mask] <- stats::ecdf(maxPostInPeak[mask])(maxPostInPeak[mask]) * 1000
    return(peakScores)
}


