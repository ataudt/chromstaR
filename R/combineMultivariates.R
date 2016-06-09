#' Combine combinatorial states from several Multivariates
#' 
#' Combine combinatorial states from several \code{\link{multiHMM}} objects. Combinatorial states can be combined for objects containing multiple marks (\code{mode='mark'}) or multiple conditions (\code{mode='condition'}).
#' 
#' @param hmms A \code{list()} with \code{\link{multiHMM}} objects. Alternatively a character vector with filenames that contain \code{\link{multiHMM}} objects.
#' @param mode Mode of combination. See \code{\link{Chromstar}} for a description of the \code{mode} parameter.
#' @return A \code{\link{combinedMultiHMM}} objects with combinatorial states for each condition.
#' @author Aaron Taudt
#' @export
#' @examples
#'### Multivariate peak calling for spontaneous hypertensive rat (SHR) ###
#'# Get example BED files for 2 different marks in hypertensive rat (SHR)
#'file.path <- system.file("extdata","euratrans", package='chromstaRData')
#'bedfiles <- list.files(file.path, full.names=TRUE, pattern='SHR')[c(1:2,6:7)]
#'# Construct experiment structure
#'exp <- data.frame(file=bedfiles, mark=c("H3K27me3","H3K27me3","H3K4me3","H3K4me3"),
#'                  condition=rep("SHR",4), replicate=c(1:2,1:2), pairedEndReads=FALSE)
#'states <- stateBrewer(exp, mode='mark')
#'# Bin the data
#'data(rn4_chrominfo)
#'binned.data <- list()
#'for (bedfile in bedfiles) {
#'  binned.data[[basename(bedfile)]] <- binReads(bedfile, binsize=1000, experiment.table=exp,
#'                                               assembly=rn4_chrominfo, chromosomes='chr12')
#'}
#'# Obtain the univariate fits
#'models <- list()
#'for (i1 in 1:length(binned.data)) {
#'  models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, eps=1)
#'}
#'# Call multivariate peaks
#'multimodel.SHR <- callPeaksMultivariate(models, use.states=states, eps=1, max.time=60)
#'
#'#'### Multivariate peak calling for brown norway (BN) rat ###
#'# Get example BED files for 2 different marks in brown norway rat
#'file.path <- system.file("extdata","euratrans", package='chromstaRData')
#'bedfiles <- list.files(file.path, full.names=TRUE, pattern='BN')[c(1:2,7:8)]
#'# Construct experiment structure
#'exp <- data.frame(file=bedfiles, mark=c("H3K27me3","H3K27me3","H3K4me3","H3K4me3"),
#'                  condition=rep("BN",4), replicate=c(1:2,1:2), pairedEndReads=FALSE)
#'states <- stateBrewer(exp, mode='mark')
#'# Bin the data
#'data(rn4_chrominfo)
#'binned.data <- list()
#'for (bedfile in bedfiles) {
#'  binned.data[[basename(bedfile)]] <- binReads(bedfile, binsize=1000, experiment.table=exp,
#'                                               assembly=rn4_chrominfo, chromosomes='chr12')
#'}
#'# Obtain the univariate fits
#'models <- list()
#'for (i1 in 1:length(binned.data)) {
#'  models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, eps=1)
#'}
#'# Call multivariate peaks
#'multimodel.BN <- callPeaksMultivariate(models, use.states=states, eps=1, max.time=60)
#'
#'### Combine multivariates ###
#'hmms <- list(multimodel.SHR, multimodel.BN)
#'comb.model <- combineMultivariates(hmms, mode='mark')
#'
combineMultivariates <- function(hmms, mode) {
    
    ## Helper function
    getCondition <- function(combinations, condition) {
        combinations <- sub('\\[','',combinations)
        combinations <- sub('\\]','',combinations)
        combinations.split <- strsplit(as.character(combinations), '\\+')
        condition.string <- paste0('\\<', condition, '\\>')
        combinations.split.condition <- lapply(combinations.split, function(x) { grep(condition.string, x, value=TRUE) })
        combinations.condition <- sapply(combinations.split.condition, paste, collapse='+')
        combinations.condition <- gsub(paste0('-',condition.string), '', combinations.condition)
        combinations.condition <- paste0('[', combinations.condition, ']')
        return(combinations.condition)
    }
      
    if (mode == 'mark') {
        ## Load first HMM for coordinates
        ptm <- startTimedMessage("Processing condition ",1," ...")
        hmm <- suppressMessages( loadHmmsFromFiles(hmms[[1]], check.class=class.multivariate.hmm)[[1]] )
        bins <- hmm$bins
        mcols(bins) <- NULL
        ## Add combinatorial states, counts and posteriors
        infos <- list()
        infos[[1]] <- hmm$info
        combs <- list()
        combs[[1]] <- hmm$bins$combination
        counts <- list()
        counts[[1]] <- hmm$bins$counts
        posteriors <- list()
        posteriors[[1]] <- hmm$bins$posteriors
        binstates <- list()
        binstates[[1]] <- dec2bin(hmm$bins$state, colnames=hmm$info$ID)
        stopTimedMessage(ptm)
        
        if (length(hmms) >= 2) {
            for (i1 in 2:length(hmms)) {
                ptm <- startTimedMessage("Processing condition ",i1," ...")
                hmm <- suppressMessages( loadHmmsFromFiles(hmms[[i1]], check.class=class.multivariate.hmm)[[1]] )
                infos[[i1]] <- hmm$info
                combs[[i1]] <- hmm$bins$combination
                counts[[i1]] <- hmm$bins$counts
                posteriors[[i1]] <- hmm$bins$posteriors
                binstates[[i1]] <- dec2bin(hmm$bins$state, colnames=hmm$info$ID)
                stopTimedMessage(ptm)
            }
        }
        counts <- do.call(cbind, counts)
        posteriors <- do.call(cbind, posteriors)
        infos <- do.call(rbind, infos)
        conditions <- unique(infos$condition)
        states <- factor(bin2dec(do.call(cbind, binstates)))
        names(states) <- NULL
        names(combs) <- conditions
        combs.df <- as(combs,'DataFrame')
        
    } else if (mode == 'condition') {
        ### Get posteriors and binary states
        infos <- list()
        counts <- list()
        posteriors <- list()
        binstates <- list()
        for (i1 in 1:length(hmms)) {
            ptm <- startTimedMessage("Processing HMM ",i1," ...")
            hmm <- suppressMessages( loadHmmsFromFiles(hmms[[i1]], check.class=class.multivariate.hmm)[[1]] )
            infos[[i1]] <- hmm$info
            counts[[i1]] <- hmm$bins$counts
            posteriors[[i1]] <- hmm$bins$posteriors
            binstates[[i1]] <- dec2bin(hmm$bins$state, colnames=hmm$info$ID)
            stopTimedMessage(ptm)
        }
        counts <- do.call(cbind, counts)
        posteriors <- do.call(cbind, posteriors)
        infos <- do.call(rbind, infos)
        conditions <- unique(infos$condition)
        binstates <- do.call(cbind, binstates)
        states <- factor(bin2dec(binstates))
        names(states) <- NULL

        bins <- hmm$bins
        mcols(bins) <- NULL

        combs <- list()
        for (condition in conditions) {
            index <- which(infos$condition==condition)
            # Make mapping
            mapping.df <- stateBrewer(infos[index,1:5], mode='mark')
            mapping <- mapping.df$combination
            names(mapping) <- mapping.df$state
            # Make combinations
            binstates.cond <- matrix(binstates[,index], ncol=length(index))
            states.cond <- factor(bin2dec(binstates.cond))
            combs[[condition]] <- mapping[as.character(states.cond)]
        }
        combs.df <- as.data.frame(combs) # get factors instead of characters
        combs.df <- as(combs.df, 'DataFrame')
        
    } else if (mode == 'full') {
        if (length(hmms) > 1) {
            stop("'hmms' must contain only one 'multiHMM' object when mode='full'.")
        }
        
        hmm <- suppressMessages( loadHmmsFromFiles(hmms[[1]], check.class=class.multivariate.hmm)[[1]] )
        bins <- hmm$bins
        mcols(bins) <- NULL
        infos <- hmm$info
        conditions <- unique(infos$condition)
        counts <- hmm$bins$counts
        posteriors <- hmm$bins$posteriors
        states <- hmm$bins$state
        combs <- list()
        for (condition in conditions) {
            ptm <- startTimedMessage("Processing condition ",condition," ...")
            mapping.condition <- getCondition(hmm$mapping, condition)
            names(mapping.condition) <- names(hmm$mapping)
          
            combs[[as.character(condition)]] <- mapping.condition[as.character(hmm$bins$state)]
            stopTimedMessage(ptm)
        }
        combs.df <- as.data.frame(combs) # get factors instead of characters
        combs.df <- as(combs.df, 'DataFrame')
        
    } else if (mode == 'replicate') {
        ## Load first HMM for coordinates
        ptm <- startTimedMessage("Processing mark-condition ",1," ...")
        hmm <- suppressMessages( loadHmmsFromFiles(hmms[[1]], check.class=class.multivariate.hmm)[[1]] )
        bins <- hmm$bins
        mcols(bins) <- NULL
        ## Add combinatorial states, counts and posteriors
        infos <- list()
        infos[[1]] <- hmm$info
        counts <- list()
        counts[[1]] <- hmm$bins$counts
        posteriors <- list()
        posteriors[[1]] <- hmm$bins$posteriors
        binstates <- list()
        binstates[[1]] <- dec2bin(hmm$bins$state, colnames=hmm$info$ID)
        stopTimedMessage(ptm)
        
        if (length(hmms) >= 2) {
            for (i1 in 2:length(hmms)) {
                ptm <- startTimedMessage("Processing mark-condition ",i1," ...")
                hmm <- suppressMessages( loadHmmsFromFiles(hmms[[i1]], check.class=class.multivariate.hmm)[[1]] )
                infos[[i1]] <- hmm$info
                counts[[i1]] <- hmm$bins$counts
                posteriors[[i1]] <- hmm$bins$posteriors
                binstates[[i1]] <- dec2bin(hmm$bins$state, colnames=hmm$info$ID)
                stopTimedMessage(ptm)
            }
        }
        counts <- do.call(cbind, counts)
        posteriors <- do.call(cbind, posteriors)
        binstates <- do.call(cbind, binstates)
        infos <- do.call(rbind, infos)
        conditions <- unique(infos$condition)
        combs <- list()
        for (condition in conditions) {
            index <- which(infos$condition==condition)
            states <- factor(bin2dec(matrix(binstates[,index], ncol=length(index))))
            names(states) <- NULL
            mapping.df <- stateBrewer(infos[index,1:5], mode='mark')
            mapping <- mapping.df$combination
            names(mapping) <- mapping.df$state
            combs[[condition]] <- mapping[as.character(states)]
        }
        names(combs) <- conditions
        combs.df <- as(combs,'DataFrame')
    } else {
        stop("Unknown mode '", mode, "'.")
    }
    # Reassign levels such that all conditions have the same levels
    comblevels <- sort(unique(unlist(lapply(combs.df, levels))))
    combs.df <- endoapply(combs.df, function(x) { x <- factor(x, levels=comblevels) })
    names(combs.df) <- paste0('combination.', names(combs.df))

    ## Assign transition groups
    ptm <- startTimedMessage("Assigning transition groups ...")
    freqs <- suppressMessages( transitionFrequencies(combstates=as.list(combs.df)) )
    bins$transition.group <- freqs$per.bin$group
    stopTimedMessage(ptm)

    ## Assign combinatorial states
    bins$state <- states
    mcols(bins)[names(combs.df)] <- combs.df
    
    ## Transferring counts and posteriors
    bins$counts <- counts
    bins$posteriors <- posteriors

    ## Add differential score ##
    bins$differential.score <- differentialScoreSum(bins$posteriors, infos)

    ### Redo the segmentation for all conditions combined
    ptm <- startTimedMessage("Redoing segmentation for all conditions combined ...")
    segments <- suppressMessages( multivariateSegmentation(bins, column2collapseBy='state') )
    names(mcols(segments)) <- setdiff(names(mcols(bins)), c('posteriors','counts'))
    stopTimedMessage(ptm)
    
    ### Redo the segmentation for each condition separately
    ptm <- startTimedMessage("Redoing segmentation for each condition separately ...")
    segments.separate <- list()
    for (cond in names(combs)) {
        bins.cond <- bins
        mcols(bins.cond) <- mcols(bins)[paste0('combination.',cond)]
        df <- as.data.frame(bins.cond)
        names(df)[6] <- cond
        segments.cond <- suppressMessages( collapseBins(df, column2collapseBy=cond, columns2drop=c('width', grep('posteriors', names(df), value=TRUE))) )
        segments.cond <- as(segments.cond, 'GRanges')
        names(mcols(segments.cond)) <- 'combination'
        seqlengths(segments.cond) <- seqlengths(bins)[seqlevels(segments.cond)]
        segments.separate[[cond]] <- segments.cond
    }
    stopTimedMessage(ptm)
    
    ### Make return object
    hmm <- list()
    class(hmm) <- class.combined.multivariate.hmm
    hmm$info <- infos
    rownames(hmm$info) <- NULL
    hmm$bins <- bins
    hmm$segments <- segments
    hmm$segments.separate <- segments.separate
    hmm$frequencies <- freqs$table
    return(hmm)
    
}
