#' Frequencies of combinatorial states
#'
#' Get the genomewide frequency of each combinatorial state.
#'
#' @param multi.hmm A \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object or a file that contains such an object.
#' @param combinations A vector with combinations for which the frequency will be calculated. If \code{NULL} all combinations will be considered.
#' @param per.mark Set to \code{TRUE} if you want frequencies per mark instead of per combination.
#' @return A table with frequencies of each combinatorial state.
#' @author Aaron Taudt
#' @export
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'genomicFrequencies(model)
#'
genomicFrequencies <- function(multi.hmm, combinations=NULL, per.mark=FALSE) {

    multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    bins <- multi.hmm$bins
    segs <- multi.hmm$segments
    peaks <- multi.hmm$peaks
    
    if (per.mark) {
        # binstates <- dec2bin(bins$state, colnames=multi.hmm$info$ID)
        # t <- colSums(binstates) / nrow(binstates)
        t <- sapply(peaks, function(peak) { sum(as.numeric(width(peak))) }) / sum(as.numeric(width(bins)))
        s <- sapply(peaks, length)
        return(list(frequency=t, domains=s))
    }
      
    if (class(multi.hmm)==class.multivariate.hmm) {

        if (is.null(combinations)) {
            comb.levels <- levels(bins$combination)
        } else {
            comb.levels <- combinations
        }
        t <- table(bins$combination) / length(bins)
        t <- t[names(t) %in% comb.levels]
        s <- table(segs$combination)
        s <- s[names(s) %in% comb.levels]
        return(list(frequency=t, domains=s))
      
    } else if (class(multi.hmm)==class.combined.multivariate.hmm) {
      
        if (is.null(combinations)) {
            comb.levels <- unique(as.vector(sapply(getCombinations(bins), levels)))
        } else {
            comb.levels <- combinations
        }
        t <- sapply(mcols(bins)[grepl('combination', names(mcols(bins)))], function(x) { table(x) / length(bins) })
        t <- t[rownames(t) %in% comb.levels,]
        s <- sapply(mcols(segs)[grepl('combination', names(mcols(segs)))], table)
        s <- s[rownames(s) %in% comb.levels,]
        return(list(frequency=t, domains=s))
      
    }
}


#' Transition frequencies of combinatorial states
#'
#' Get a table of transition frequencies between combinatorial states of different \code{\link{multiHMM}}s.
#'
#' @param multi.hmms A named list with \code{\link{multiHMM}} objects or a vector with filenames that contain such objects.
#' @param combined.hmm A \code{\link{combinedMultiHMM}} object. If specified, \code{multi.hmms} is ignored.
#' @param zero.states The string(s) which identifies the zero.states.
#' @param combstates Alternative input instead of \code{multi.hmms}: A named list of combinatorial state vectors instead of HMMs. Names must be of the form "combination.X", where X is an arbitrary string. If this is specified, \code{multi.hmms} and \code{combined.hmm} will be ignored.
#' @return A data.frame with transition frequencies.
#' @author Aaron Taudt
#' @export
#' @examples 
#'#=== Step 1: Preparation ===
#'## Prepare the file paths. Exchange this with your input and output directories.
#'inputfolder <- system.file("extdata","euratrans", package="chromstaRData")
#'outputfolder <- file.path(tempdir(), 'SHR-BN-example')
#'## Define experiment structure
#'data(experiment_table)
#'print(experiment_table)
#'## Define assembly
#'# This is only necessary if you have BED files, BAM files are handled automatically.
#'# For common assemblies you can also specify them as 'hg19' for example.
#'data(rn4_chrominfo)
#'head(rn4_chrominfo)
#'
#'#=== Step 2: Run Chromstar ===
#'## Run ChromstaR
#'Chromstar(inputfolder, experiment.table=experiment_table,
#'          outputfolder=outputfolder, numCPU=2, binsize=1000, assembly=rn4_chrominfo,
#'          prefit.on.chr='chr12', chromosomes='chr12', mode='combinatorial', eps.univariate=1,
#'          eps.multivariate=1)
#'## Results are stored in 'outputfolder' and can be loaded for further processing
#'list.files(outputfolder)
#'model <- get(load(file.path(outputfolder,'combined', 'combined_mode-combinatorial.RData')))
#'
#'#=== Step 3: Analysis ===
#'# Get frequencies
#'freqs <- transitionFrequencies(combined.hmm=model)
#'freqs$table
#'
transitionFrequencies <- function(multi.hmms=NULL, combined.hmm=NULL, zero.states="[]", combstates=NULL) {

    if (is.null(combstates)) {
        if (is.null(combined.hmm)) {
            if (is.null(names(multi.hmms))) {
                stop("'multi.hmms' must be a named list of multiHMM objects.")
            }
            if (!all(grepl('^combination.', names(multi.hmms)))) {
                stop("'multi.hmms' must be a named list of multiHMM objects. Names must have the form 'combination.X', where X is an arbitrary string.")
            }
            ## Get combinatorial states in loop to save memory
            ptm <- startTimedMessage("Loading HMMs ...")
            combstates <- list()
            for (imodel in 1:length(multi.hmms)) {
                multi.hmm <- suppressMessages( loadHmmsFromFiles(multi.hmms[[imodel]], check.class=class.multivariate.hmm)[[1]] )
                combstates[[imodel]] <- multi.hmm$bins$combination
            }
            names(combstates) <- names(multi.hmms)
            stopTimedMessage(ptm)
        } else {
            combined.hmm <- suppressMessages( loadHmmsFromFiles(combined.hmm, check.class=class.combined.multivariate.hmm)[[1]] )
            combstates <- as.list(mcols(combined.hmm$bins)[grepl('combination', names(mcols(combined.hmm$bins)))])
        }
    } else {
        if (is.null(names(combstates))) {
            stop("'combstates' must be a named list.")
        }
        if (!all(grepl('^combination.', names(combstates)))) {
            stop("'combstates' must be a named list. Names must have the form 'combination.X', where X is an arbitrary string.")
        }
    }
    num.models <- length(combstates)
    conditions <- names(combstates)

    ### Get transitions for whole genome ###
    ptm <- startTimedMessage("Getting transitions ...")
    combstates$sep <- '<>'
    gentrans <- do.call(paste, combstates)
    tab <- table(gentrans) / length(gentrans)
    trans <- matrix(unlist(strsplit(names(tab), '<>')), ncol=length(conditions), dimnames=list(transition=NULL, condition=conditions), byrow = TRUE)
    freqtrans <- data.frame(trans)
    freqtrans$frequency <- as.numeric(tab)
    freqtrans$transition <- names(tab)
    freqtrans <- freqtrans[order(freqtrans$frequency, decreasing=TRUE),]
    rownames(freqtrans) <- NULL
    # Cumulative frequencies
    freqtrans$cumulative.frequency <- cumsum(freqtrans$frequency)
    stopTimedMessage(ptm)
    
    ### Number of domains ###
    ptm <- startTimedMessage("Number of domains ...")
    rle.gentrans <- rle(as.character(gentrans))
    ndomains <- table(rle.gentrans$values)
    freqtrans$domains <- ndomains[freqtrans$transition]
    stopTimedMessage(ptm)

    ### Assigning groups for frequency table ###
    ptm <- startTimedMessage("Assigning groups ...")
    freqtrans <- assignGroups(freqtrans, zero.states, num.models)
    stopTimedMessage(ptm)

    ### Assigning groups over whole genome ###
    ptm <- startTimedMessage("Assigning groups for whole genome ...")
    mapping <- freqtrans$group
    names(mapping) <- freqtrans$transition
    gengroups <- mapping[gentrans]
    gentrans <- data.frame(transition=gentrans, group=gengroups)
    stopTimedMessage(ptm)

    ## Remove unneeded column
    freqtrans$transition <- NULL
    
    ## Reorder columns
    freqtrans <- freqtrans[, c(grep('combination', names(freqtrans), value=TRUE), 'domains', 'frequency', 'cumulative.frequency', 'group')]
    
    ## Return value ##
    return(list(table=freqtrans, per.bin=gentrans))

}


assignGroups <- function(freqtrans, zero.states, num.models) {
    
    freqtrans$group <- 'other'
    # Stage-specific and constant states
    levels.combstates <- unique(unlist(lapply(freqtrans[,grep('combination', names(freqtrans))], levels)))
    levels.combstates <- setdiff(levels.combstates, zero.states)
    levels.combstates <- gsub('\\+','\\\\+',levels.combstates)
    levels.combstates <- gsub('\\[','\\\\[', levels.combstates)
    levels.combstates <- gsub('\\]','\\\\]', levels.combstates)
    ## Determine number of zero.states in each transition
    combs <- as.matrix(freqtrans[,grep("combination", names(freqtrans))])
    zero.matrix <- apply(combs, 2, function(x) { x %in% zero.states })
    num.zeros <- rowSums(zero.matrix)
    ## Stage specific transitions
    mask <- num.zeros == (num.models - 1)
    combination <- t(combs)[,mask][!t(zero.matrix)[,mask]]
    freqtrans$group[mask] <- paste0('stage-specific ', gsub('\\\\','',combination))
    ## Constant transitions
    mask <- rep(TRUE, nrow(freqtrans))
    if (ncol(combs) >= 2) {
        for (i1 in 2:ncol(combs)) {
            mask <- mask & (combs[,1] == combs[,i1])
        }
    }
    combination <- combs[mask,1]
    freqtrans$group[mask] <- paste0('constant ', gsub('\\\\','',combination))
    # Zero transitions
    mask <- num.zeros == num.models
    freqtrans$group[mask] <- 'zero transition'
    return(freqtrans)
}
