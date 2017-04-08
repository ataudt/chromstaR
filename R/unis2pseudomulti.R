#' Combine univariate HMMs to a multivariate HMM
#'
#' Combine multiple \code{\link{uniHMM}}s to a \code{\link{multiHMM}} without running \code{\link{callPeaksMultivariate}}. This should only be done for comparison purposes.
#'
#' Use this function if you want to combine ChIP-seq samples without actually running a multivariate Hidden Markov Model. The resulting object will be of class \code{\link{multiHMM}} but will not be truly multivariate.
#'
#' @author Aaron Taudt
#' @param hmms A named list of \code{\link{uniHMM}} objects. Names will be used to generate the combinations.
#' @return A \code{\link{multiHMM}} object.
#' @export
#' @examples
#'# Get example BAM files for 2 different marks in hypertensive rat (SHR)
#'file.path <- system.file("extdata","euratrans", package='chromstaRData')
#'files <- list.files(file.path, full.names=TRUE, pattern='SHR.*bam$')[c(1,4)]
#'# Bin the data
#'data(rn4_chrominfo)
#'binned.data <- list()
#'for (file in files) {
#'  binned.data[[basename(file)]] <- binReads(file, binsizes=1000, stepsizes=500,
#'                                            assembly=rn4_chrominfo, chromosomes='chr12')
#'}
#'# Obtain the univariate fits
#'models <- list()
#'for (i1 in 1:length(binned.data)) {
#'  models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, eps=1)
#'}
#'## Combine the univariate HMMs without fitting a multivariate HMM
#'names(models) <- c('H3K27me3','H3K4me3')
#'pseudo.multi.HMM <- unis2pseudomulti(models)
#'## Compare frequencies with real multivariate HMM
#'exp <- data.frame(file=files, mark=c("H3K27me3","H3K4me3"),
#'                  condition=rep("SHR",2), replicate=c(1,1), pairedEndReads=FALSE,
#'                  controlFiles=NA)
#'states <- stateBrewer(exp, mode='combinatorial')
#'real.multi.HMM <- callPeaksMultivariate(models, use.states=states, eps=1, max.time=60)
#'genomicFrequencies(real.multi.HMM)
#'genomicFrequencies(pseudo.multi.HMM)
#'
unis2pseudomulti <- function(hmms) {

    # Load models
    hmms <- loadHmmsFromFiles(hmms, check.class=class.univariate.hmm)

    # Extract coordinates and other stuff
    nummod = length(hmms)
    bins <- hmms[[1]]$bins
    bins$counts.rpkm <- NULL
    bins$state <- NULL
    numbins = length(hmms[[1]]$bins)
    info <- do.call(rbind, lapply(hmms, function(x) { x$info }))
    distributions = lapply(hmms,"[[","distributions")
    weights = lapply(hmms,"[[","weights")

    # Extract the reads
    ptm <- startTimedMessage("Extracting read counts from hmms ...")
    reads = matrix(NA, ncol=nummod, nrow=numbins)
    colnames(reads) <- info$ID
    for (imod in 1:nummod) {
        reads[,imod] = hmms[[imod]]$bins$counts.rpkm
    }
    bins$counts.rpkm <- reads
    stopTimedMessage(ptm)

    ## Get combinatorial states
    ptm <- startTimedMessage("Getting combinatorial states ...")
    combstates.per.bin = combinatorialStates(hmms)
    comb.states.table = table(combstates.per.bin)
    comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
    numstates <- length(comb.states)
    bins$state <- factor(combstates.per.bin, levels=comb.states)
    binary.comb.states <- dec2bin(comb.states, colnames=names(hmms))
    binary.comb.states.list <- list()
    for (icol in 1:ncol(binary.comb.states)) {
        binary.comb.states.list[[colnames(binary.comb.states)[icol]]] <- c('',colnames(binary.comb.states)[icol])[binary.comb.states[,icol]+1]
    }
    binary.comb.states.list$sep='+'
    mapping <- do.call(paste, binary.comb.states.list)
    mapping <- gsub('\\+{2,}', '+', mapping)
    mapping <- sub('^\\+', '', mapping)
    mapping <- sub('\\+$', '', mapping)
    mapping <- paste0('[', mapping, ']')
    names(mapping) <- rownames(binary.comb.states)
    bins$combination <- factor(mapping[as.character(bins$state)], levels=mapping)
    stopTimedMessage(ptm)
    
    ## Calculate transition matrix
    ptm <- startTimedMessage("Estimating transition matrix ...")
    df <- data.frame(from = combstates.per.bin[-length(combstates.per.bin)], to = combstates.per.bin[-1])
    t <- table(df)
    t <- sweep(t, 1, rowSums(t), "/")
    # Select only states that are in data
    t <- t[as.character(comb.states),as.character(comb.states)]
    t.estimated <- t
    stopTimedMessage(ptm)

    ## Return multi.hmm
    result <- list()
    result$info <- info
    rownames(result$info) <- NULL
    result$bins <- bins
    ## Segmentation
        df <- as.data.frame(result$bins)
        ind.readcols <- grep('^counts', names(df))
        ind.widthcol <- grep('width', names(df))
        ind.scorecol <- grep('score', names(df))
        red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c(ind.readcols, ind.widthcol)))
        red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'], combination=red.df[,'combination'])
        result$segments <- red.gr
        seqlengths(result$segments) <- seqlengths(result$bins)[seqlevels(result$segments)]
    ## Parameters
        result$mapping <- mapping
        # Weights
        tstates <- table(combstates.per.bin)
        result$weights <- sort(tstates/sum(tstates), decreasing=TRUE)
        # Transition matrices
        result$transitionProbs <- t.estimated
        # Distributions
        result$distributions <- distributions
        names(result$distributions) <- info$ID
    ## Convergence info
        convergenceInfo <- list(eps=Inf, loglik=Inf, loglik.delta=Inf, num.iterations=Inf, time.sec=Inf)
        result$convergenceInfo <- convergenceInfo
    ## Correlation matrices
#         result$correlation.matrix <- correlationMatrix2use
    ## Add class
        class(result) <- class.multivariate.hmm

    return(result)

}
