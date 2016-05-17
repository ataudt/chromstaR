#' Frequencies of combinatorial states
#'
#' Get the genomewide frequency of each combinatorial state.
#'
#' @param multi.hmm A \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object or a file that contains such an object.
#' @param combinations A vector with combinations for which the frequency will be calculated. If \code{NULL} all combinations will be considered.
#' @return A table with frequencies of each combinatorial state.
#' @author Aaron Taudt
#' @export
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-mark_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'genomicFrequencies(model)
#'
genomicFrequencies <- function(multi.hmm, combinations=NULL) {

    multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    bins <- multi.hmm$bins
      
    if (class(multi.hmm)==class.multivariate.hmm) {

        if (is.null(combinations)) {
            comb.levels <- levels(bins$combination)
        } else {
            comb.levels <- combinations
        }
        t <- table(bins$combination) / length(bins)
        t <- t[names(t) %in% comb.levels]
        return(t)
      
    } else if (class(multi.hmm)==class.combined.multivariate.hmm) {
      
        if (is.null(combinations)) {
            comb.levels <- unique(as.vector(sapply(mcols(bins), levels)))
        } else {
            comb.levels <- combinations
        }
        t <- sapply(mcols(bins), function(x) { table(x) / length(bins) })
        t <- t[rownames(t) %in% comb.levels,]
        return(t)
      
    }
}


#' Transition frequencies of combinatorial states
#'
#' Get a table of transition frequencies between combinatorial states of different \code{\link{multiHMM}}s.
#'
#' @param multi.hmms A named list with \code{\link{multiHMM}} objects or a vector with filenames that contain such objects.
#' @param combined.hmm A \code{\link{combinedMultiHMM}} object. If specified, \code{multi.hmms} is ignored.
#' @param zero.states The string(s) which identifies the zero.states.
#' @param combstates Alternative input instead of \code{multi.hmms}: A named list of combinatorial state vectors instead of HMMs. If this is specified, \code{multi.hmms} and \code{combined.hmm} will be ignored.
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

#'#=== Step 2: Run Chromstar ===
#'## Run ChromstaR
#'Chromstar(inputfolder, experiment.table=experiment_table,
#'          outputfolder=outputfolder, numCPU=2, binsize=1000, assembly=rn4_chrominfo,
#'          prefit.on.chr='chr12', mode='mark', eps.univariate=1, eps.multivariate=1)
#'## Results are stored in 'outputfolder' and can be loaded for further processing
#'list.files(outputfolder)
#'model <- get(load(file.path(outputfolder,'combined', 'combined_mode-mark.RData')))

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
    for (combination in levels.combstates) {
        string.other.levels <- paste(setdiff(levels.combstates,combination), collapse='|')
        if (string.other.levels=="") {
            string.other.levels <- 'AAAAAA' # workaround
        }
        mask <- intersect(grep(combination, freqtrans$transition), grep(string.other.levels, freqtrans$transition, invert=TRUE))
        freqtrans$group[mask] <- paste0('stage-specific ',gsub('\\\\','',combination))
        mask <- sapply(gregexpr(combination,freqtrans$transition), function(x) { length(which(x!=-1)) }) == num.models
        freqtrans$group[mask] <- paste0('constant ', gsub('\\\\','',combination))
    }
    # Zero transitions
    freqtrans.split <- strsplit(sub('<>$','<><>',as.character(freqtrans$transition)),'<>')
    freqtrans.split <- do.call(rbind, freqtrans.split)
    df <- as.data.frame(apply(freqtrans.split, 2, function(x) { x %in% zero.states } ))
    # if (ncol(df)==1) {
    #     iszero <- Reduce('&', df[,1])
    # } else {
        iszero <- Reduce('&', as.list(df))
    # }
    freqtrans$group[iszero] <- 'zero transition'
    return(freqtrans)
}
