#' Simulate univariate data
#' 
#' Simulate known states, read counts and read coordinates using a univariate Hidden Markov Model with three states ("zero-inflation", "unmodified" and "modified").
#' 
#' @param bins A \code{\link[GenomicRanges]{GRanges}} object for which reads will be simulated.
#' @param transition A matrix with transition probabilities.
#' @param emission A data.frame with emission distributions (see \code{\link{uniHMM}} entry 'distributions').
#' @inheritParams simulateReadsFromCounts
#' @return A \code{list} with entries $bins containing the simulated states and read count, $reads with simulated read coordinates and $transition and $emission.
#' @importFrom stats runif rnbinom aggregate
simulateUnivariate <- function(bins, transition, emission, fragLen=50) {

    # Calculate some variables
    numstates = ncol(transition)
    if (numstates!=3) {
        stop("The transition matrix is expected to have 3 columns and 3 rows")
    }
    numbins = length(bins)

    ## Make state vector from transition matrix
    # Generate sample of random numbers
    rsample = stats::runif(numbins,0,1)
    # Integrate transition matrix by row and add -1 in front
    cumtransition = cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
    # Generate the state vector by going through each state
    ptm <- startTimedMessage("Generating states from transition matrix ...")
    states = matrix(rep(NA,numstates*numbins), ncol=numstates)
    for (irow in 1:numstates) {
        states[,irow] = findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
    }
    statevec = rep(NA,numbins)
    statevec[1] = 1
    for (ibin in 2:numbins) {
        statevec[ibin] = states[ibin,statevec[ibin-1]]
    }
    stopTimedMessage(ptm)

    ## Make reads from state vector and distributions
    # Generate the read counts by drawing from distribution
    ptm <- startTimedMessage("Generating read counts from emission parameters and states ...")
    reads = rep(NA, numbins)
    numbins.in.state = stats::aggregate(rep(1,length(statevec)), list(state=statevec), sum)
    reads[statevec==1] = 0
    if (!is.na(numbins.in.state[2,'x'])) {
        reads[statevec==2] = stats::rnbinom(numbins.in.state[2,'x'], size=emission[2,'size'], prob=emission[2,'prob'])
    }
    if (!is.na(numbins.in.state[3,'x'])) {
        reads[statevec==3] = stats::rnbinom(numbins.in.state[3,'x'], size=emission[3,'size'], prob=emission[3,'prob'])
    }
    stopTimedMessage(ptm)

    ## Combine in GRanges
    bins.out <- bins
    mcols(bins.out) <- NULL
    bins.out$counts <- reads
    bins.out$state <- state.labels[statevec]

    # Return the output
    out = list(bins = bins.out,
               reads = simulateReadsFromCounts(bins.out),
               transition = transition,
               emission = emission
               )
    return(out)

}


#' Simulate read coordinates
#' 
#' Simulate read coordinates using read counts as input.
#' 
#' @param bins A \code{\link[GenomicRanges]{GRanges}} with read counts.
#' @param fragLen Length of the simulated read fragments.
#' @return A \code{\link[GenomicRanges]{GRanges}} with read coordinates.
simulateReadsFromCounts <- function(bins, fragLen=50) {
    
    ptm <- startTimedMessage("Generating read coordinates ...")
    ## Do this for every chromosome
    bins.splt <- split(bins, seqnames(bins))
    reads <- GRangesList()
    seqlevels(reads) <- seqlevels(bins)
    seqlengths(reads) <- seqlengths(bins)
    for (chrom in names(bins.splt)) {
        bins.chrom <- bins.splt[[chrom]]
        if (length(bins.chrom) > 0) {
            total.reads <- sum(bins.chrom$counts)
            binsize <- mean(width(bins.chrom))
            ## For each read, get the start position
            rlcounts <- rle(1)
            rlcounts$lengths <- bins.chrom$counts
            rlcounts$values <- as.numeric(start(bins.chrom))
            starts <- inverse.rle(rlcounts)
            ## For each read, simulate a random offset to the start of the bin
            offsets <- round(runif(total.reads, min=0, max=binsize))
            strand <- c('+','-')[round(runif(total.reads, min=1, max=2))]
            ## Make read coordinates
            if (length(starts) > 0) {
                reads[[chrom]] <- sort(GRanges(seqnames=chrom, ranges=IRanges(start=starts+offsets, end=starts+offsets+fragLen)))
            } else {
                reads[[chrom]] <- GRanges()
            }
            strand(reads[[chrom]]) <- strand
        }
    }
    reads <- unlist(reads, use.names=FALSE)
    stopTimedMessage(ptm)
    
    return(reads)
    
}


#' Simulate multivariate data
#' 
#' Simulate known states, read counts and read coordinates using a multivariate Hidden Markov Model.
#' 
#' @param bins A \code{\link[GenomicRanges]{GRanges}} object for which reads will be simulated.
#' @param transition A matrix with transition probabilities.
#' @param emissions A list() with data.frames with emission distributions (see \code{\link{uniHMM}} entry 'distributions').
#' @param weights A list() with weights for the three univariate states.
#' @param correlationMatrices A list with correlation matrices.
#' @param combstates A vector with combinatorial states.
#' @param IDs A character vector with IDs.
#' @inheritParams simulateReadsFromCounts
#' @return A \code{list()} with entries $bins containing the simulated states and read count, $reads with simulated read coordinates.
#' @importFrom mvtnorm rmvnorm
# #' @examples 
# #'#'### Get an example multiHMM ###
# #'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
# #'                     package="chromstaR")
# #'model <- get(load(file))
# #'### Simulate reads ###
# #'sim <- chromstaR:::simulateMultivariate(bins=model$bins,
# #'                     transition=model$transitionProbs,
# #'                     emission=model$distributions,
# #'                     weights=model$weights.univariate,
# #'                     correlationMatrices=model$correlation.matrix,
# #'                     combstates=names(model$mapping),
# #'                     IDs=model$info$ID)
simulateMultivariate = function(bins, transition, emissions, weights, correlationMatrices, combstates, IDs, fragLen=50) {

    # Calculate some variables
    numstates <- ncol(transition)
    numbins <- length(bins)
    numtracks <- length(emissions)

    ## Make state vector from transition matrix
    # Generate sample of random numbers
    rsample <- stats::runif(numbins,0,1)
    # Integrate transition matrix by row and add -1 in front
    cumtransition <- cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
    # Generate the state vector by going through each state
    ptm <- startTimedMessage("Generating states from transition matrix...")
    states <- matrix(NA, ncol=numstates, nrow=numbins)
    for (irow in 1:numstates) {
        states[,irow] <- findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
    }
    statevec <- rep(NA,numbins)
    statevec[1] <- 1
    for (ibin in 2:numbins) {
        statevec[ibin] <- states[ibin,statevec[ibin-1]]
    }
    # Replace the states by combinatorial states
    statevec <- combstates[statevec]
    stopTimedMessage(ptm)

    ## Make reads from state vector and emission distributions
    counts <- matrix(rep(NA, numtracks*numbins), ncol=numtracks)
    colnames(counts) <- IDs

    sizes <- unlist(lapply(emissions, "[", 2:3, 'size'))
    probs <- unlist(lapply(emissions, "[", 2:3, 'prob'))
    ws1 <- unlist(lapply(weights,"[",1))
    ws2 <- unlist(lapply(weights,"[",2))
    ws3 <- unlist(lapply(weights,"[",3))
    ws <- ws1 / (ws2+ws1)
    for (istate in combstates) {
        message("Generating counts for state ",istate)
        i <- which(combstates==istate)

        # Convert istate to binary representation
        binary_state <- rev(as.integer(intToBits(istate))[1:numtracks])
        binary_state <- dec2bin(istate, colnames=IDs)
        binary_stateTF <- unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])

        # Draw from the multivariate normal
        n <- length(which(statevec==istate))
        if (n > 0) {
            ptm <- startTimedMessage("  Drawing from multivariate normal ...")
            z <- matrix(mvtnorm::rmvnorm(n, mean=rep(0,numtracks), sigma=correlationMatrices[,,i]), ncol=numtracks)
            stopTimedMessage(ptm)
            # Transform to uniform space
            ptm <- startTimedMessage("  Transforming to uniform space ...")
            u <- matrix(apply(z, 2, stats::pnorm), ncol=numtracks)
            stopTimedMessage(ptm)
            # Transform to count space using marginals
            ptm <- startTimedMessage("  Transforming to count space ...")
            isizes <- sizes[binary_stateTF]
            iprobs <- probs[binary_stateTF]
            for (imod in 1:numtracks) {
                mask <- statevec==istate
                if (binary_state[imod] == 0) {
                    counts[mask,imod] <- qzinbinom(u[,imod], w=ws[imod], size=isizes[imod], prob=iprobs[imod])
                } else {
                    counts[mask,imod] <- stats::qnbinom(u[,imod], size=isizes[imod], prob=iprobs[imod])
                }
            }
            stopTimedMessage(ptm)
        }

    }

    ## Combine in GRanges
    bins.out <- bins
    mcols(bins.out) <- NULL
    bins.out$counts <- counts
    bins.out$state <- factor(statevec, levels=combstates)

    ## Make reads
    reads <- list()
    for (icol in 1:ncol(bins.out$counts)) {
        track <- colnames(bins.out$counts)[icol]
        bins.temp <- bins.out
        bins.temp$counts <- bins.out$counts[,icol]
        reads[[track]] <- simulateReadsFromCounts(bins.temp)
    }
    
    # Return the output
    out <- list(bins = bins.out,
                reads = reads,
                emissions = emissions,
                transition = transition,
                correlationMatrices = correlationMatrices,
                combstates = combstates,
                weights = weights,
                IDs = IDs
                )
    return(out)

}

