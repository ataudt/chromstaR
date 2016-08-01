#' Obtain combinatorial states from experiment table
#'
#' This function computes combinatorial states from an \code{\link{experiment.table}}.
#'
#' The binary modification state (unmodified=0 or modified=1) of multiple ChIP-seq samples defines a (decimal) combinatorial state such as:
#' \tabular{ccccccc}{
#'  \tab sample1 \tab sample2 \tab sample3 \tab sample4 \tab sample5 \tab combinatorial state \cr
#' bin1 \tab 0 \tab 0 \tab 1 \tab 0 \tab 0 \tab 4 \cr
#' bin2 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
#' bin3 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \tab 10 \cr
#' bin4 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 15 \cr
#' bin5 \tab 0 \tab 0 \tab 1 \tab 0 \tab 1 \tab 5 \cr
#' }
#'
#' @author Aaron Taudt
#' @param experiment.table A \code{data.frame} specifying the experiment structure. See \code{\link{experiment.table}}.
#' @inheritParams state.brewer
#' @param mode Mode of brewing. See \code{\link{Chromstar}} for a description of the parameter.
#' @return A data.frame with combinations and their corresponding (decimal) combinatorial states.
#' @export
#' @examples 
#'## Construct an experiment table
#'data(experiment_table)
#'print(experiment_table)
#'## Construct combinatorial states
#'stateBrewer(experiment_table, mode='combinatorial')
#'stateBrewer(experiment_table, mode='differential')
#'stateBrewer(experiment_table, mode='full', common.states=TRUE)
#'
#'## Exclude states with exclusive.table
#'excl <- data.frame(mark=c('H3K4me3','H3K27me3'),
#'                              group=c(1,1))
#'stateBrewer(experiment_table, mode='full', exclusive.table=excl)
#'
stateBrewer <- function(experiment.table, mode, differential.states=FALSE, common.states=FALSE, exclusive.table=NULL, binary.matrix=NULL) {

    check.experiment.table(experiment.table)
    exp <- experiment.table
    if (mode == 'full') {
        combstates <- state.brewer(replicates=paste0(exp$mark, '-', exp$condition), conditions=exp$condition, tracks2compare=exp$mark, differential.states=differential.states, common.states=common.states, exclusive.table=exclusive.table, binary.matrix=binary.matrix)
    } else if (mode == 'combinatorial') {
        combstates <- state.brewer(replicates=exp$mark, conditions=exp$condition, tracks2compare=exp$mark, differential.states=differential.states, common.states=common.states, exclusive.table=exclusive.table, binary.matrix=binary.matrix)
    } else if (mode == 'differential') {
        combstates <- state.brewer(replicates=exp$condition, conditions=exp$condition, tracks2compare=exp$mark, differential.states=differential.states, common.states=common.states, exclusive.table=exclusive.table, binary.matrix=binary.matrix)
    } else {
        stop("Unknown mode.")
    }
        
    return(combstates)
}


#' Obtain combinatorial states from specification
#'
#' This function returns all combinatorial (decimal) states that are consistent with a given abstract specification.
#'
#' The binary modification state (unmodified=0 or modified=1) of multiple ChIP-seq samples defines a (decimal) combinatorial state such as:
#' \tabular{ccccccc}{
#'  \tab sample1 \tab sample2 \tab sample3 \tab sample4 \tab sample5 \tab combinatorial state \cr
#' bin1 \tab 0 \tab 0 \tab 1 \tab 0 \tab 0 \tab 4 \cr
#' bin2 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
#' bin3 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \tab 10 \cr
#' bin4 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 15 \cr
#' bin5 \tab 0 \tab 0 \tab 1 \tab 0 \tab 1 \tab 5 \cr
#' }
#'
#' @author Aaron Taudt, David Widmann
#' @param replicates A vector specifying the replicate structure. Similar entries will be treated as replicates.
#' @param differential.states A logical specifying whether differential states shall be returned.
#' @param min.diff The minimum number of differences between conditions.
#' @param common.states A logical specifying whether common states shall be returned.
#' @param conditions A vector with the same length as \code{replicates}. Similar entries will be treated as belonging to the same condition. Usually your tissue or cell types or time points.
#' @param tracks2compare A vector with the same length as \code{replicates}. This vector defines the tracks between which conditions are compared. Usually your histone marks.
#' @param sep Separator used to separate the tracknames in the combinations. The default '+' should not be changed because it is assumed in follow-up functions.
#' @param statespec If this parameter is specified, \code{replicates} will be ignored. A vector composed of any combination of the following entries: \code{'0.[]', '1.[]', 'x.[]', 'r.[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'0.A'}: sample A is 'unmodified'
#'     \item \code{'1.B'}: sample B is 'modified'
#'     \item \code{'x.C'}: sample C can be both 'unmodified' or 'modified'
#'     \item \code{'r.D'}: all samples in group D have to be in the same state
#'     \item \code{'r.[]'}: all samples in group [] have to be in the same state
#'   }
#' @param diffstatespec A vector composed of any combination of the following entries: \code{'x.[]', 'd.[]'}, where [] can be any string.
#'   \itemize{
#'     \item \code{'x.A'}: sample A can be both 'unmodified' or 'modified'
#'     \item \code{'d.B'}: at least one sample in group B has to be different from the other samples in group A 
#'     \item \code{'d[]'}: at least one sample in group [] has to be different from the other samples in group [] 
#'   }
#' 
#' @param exclusive.table A \code{data.frame} or tab-separated text file with columns 'mark' and 'group'. Histone marks with the same group will be treated as mutually exclusive.
#' @param binary.matrix A logical matrix produced by \code{\link{dec2bin}}. If this is specified, only states specified by the rows of this matrix will be considered. The number of columns must match \code{length(replicates)} or \code{length(statespec)}. Only for advanced use. No error handling for incorrect input.
#' @return A data.frame with combinations and their corresponding (decimal) combinatorial states.
#' @examples
#'# Get all combinatorial states where sample1=0, sample2=1, sample3=(0 or 1),
#'#  sample4=sample5
#'chromstaR:::state.brewer(statespec=c('0.A','1.B','x.C','r.D','r.D'))
#'
#'# Get all combinatorial states where sample1=sample2=sample3, sample4=sample5
#'chromstaR:::state.brewer(statespec=c('r.A','r.A','r.A','r.B','r.B'))
#'
#'# Get all combinatorial states where sample1=sample5, sample2=sample3=1,
#'#  sample4=(0 or 1)
#'chromstaR:::state.brewer(statespec=c('r.A','1.B','1.C','x.D','r.A'))
#'
state.brewer <- function(replicates=NULL, differential.states=FALSE, min.diff=1, common.states=FALSE, conditions=NULL, tracks2compare=NULL, sep='+', statespec=NULL, diffstatespec=NULL, exclusive.table=NULL, binary.matrix=NULL) {

#     ## Debug
# #     conditions <- tissues
#     conditions <- strains
#     tracks2compare <- marks
#     differential.states <- TRUE
#     common.states <- FALSE
#     min.diff <- 1
#     statespec <- NULL
#     diffstatespec <- NULL
#     sep='+'
#     replicates <- c("Bre.H3K27Ac", "Bre.H3K27me3", "Bre.H3K4me3", "Bre.H4K20me1", "Bre.H3K27Ac", "Bre.H3K27me3", "Bre.H3K4me3", "Bre.H4K20me1", "Bre.H3K27Ac", "Bre.H3K27me3", "Bre.H3K4me3", "Bre.H4K20me1", "Gua.H3K27Ac", "Gua.H3K27me3", "Gua.H3K4me3", "Gua.H4K20me1", "Gua.H3K27Ac", "Gua.H3K27Ac", "Gua.H3K27me3", "Gua.H3K4me3", "Gua.H4K20me1")
#     conditions <- c("Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Bre", "Gua", "Gua", "Gua", "Gua", "Gua", "Gua", "Gua", "Gua", "Gua")
#     tracks2compare <- c("H3K27Ac", "H3K27me3", "H3K4me3", "H4K20me1", "H3K27Ac", "H3K27me3", "H3K4me3", "H4K20me1", "H3K27Ac", "H3K27me3", "H3K4me3", "H4K20me1", "H3K27Ac", "H3K27me3", "H3K4me3", "H4K20me1", "H3K27Ac", "H3K27Ac", "H3K27me3", "H3K4me3", "H4K20me1")
#     statespec <- paste0('r.', replicates)

    ## Check user input
    if (is.null(statespec)) {
        if (!is.null(replicates)) {
            statespec <- paste0('r.', replicates)
        } else {
            stop("Please specify either 'replicates' or 'statespec'.")
        }
    }
    for (spec in statespec) {
        if (!grepl('^1\\.', spec) & !grepl('^0\\.', spec) & !grepl('^x\\.', spec) & !grepl('^r\\.', spec)) {
            stop("argument 'statespec' expects a vector composed of any combination of the following entries: '1.[]','0.[]','x.[]','r.[]', where [] can be any string.")
        }
    }
    if (!is.null(diffstatespec)) {
        for (spec in diffstatespec) {
            if (!grepl('^1\\.', spec) & !grepl('^0\\.', spec) & !grepl('^x\\.', spec) & !grepl('^d\\.', spec)) {
                stop("argument 'diffstatespec' expects a vector composed of any combination of the following entries: '1.[]','0.[]','x.[]','r.[]', where [] can be any string.")
            }
        }
        if (length(statespec)!=length(diffstatespec)) {
            stop("argument 'diffstatespec' must have the same number of elements as 'statespec'")
        }
    }
    if (differential.states | common.states) {
        if (is.null(conditions) | is.null(tracks2compare)) {
            stop("Please specify 'conditions' and 'tracks2compare' if you want to obtain differential or common states.")
        }
    }
    if (!is.null(binary.matrix)) {
        if (class(binary.matrix) != 'matrix' | mode(binary.matrix) != 'logical') {
            stop("argument 'binary.matrix' expects a logical matrix")
        }
    }

    ## Variables
    tracknames <- sub('^.\\.', '', statespec)

    ### Generate replicate-reduced binary states ###
    numtracks <- length(statespec)
    groups <- levels(factor(statespec))
    
    if (is.null(binary.matrix)) {
        numstates <- 2^(length(which(grepl('^r\\.', groups))) + length(which(grepl('^x\\.', statespec))))
        binstates <- matrix(FALSE, ncol=numtracks, nrow=numstates)
        i1 <- 1
        for (group in groups) {
            track.index <- which(statespec==group)
            for (itrack in track.index) {
                if (grepl('^1\\.', group)) {
                    binstates[,itrack] <- TRUE
                } else if (grepl('^x\\.', group)) {
                    numeach <- numstates/2 / 2^(i1-1)
                    binstates[,itrack] <- rep(c(rep(FALSE, numeach), rep(TRUE, numeach)), 2^(i1-1))
                    i1 <- i1 + 1
                } else if (grepl('^r\\.', group)) {
                    numeach <- numstates/2 / 2^(i1-1)
                    binstates[,itrack] <- rep(c(rep(FALSE, numeach), rep(TRUE, numeach)), 2^(i1-1))
                }
            }
            if (grepl('^r\\.', group)) {
                i1 <- i1 + 1
            }
        }
    } else {
        if (ncol(binary.matrix) != length(tracknames)) {
            stop("The number of columns in 'binary.matrix' must match the number of 'replicates' or 'statespec'.")
        }
        binstates <- binary.matrix
    }
    colnames(binstates) <- tracknames

    ## Construct state names
    mask <- !grepl('^r\\.', statespec) | !duplicated(statespec) # only first replicate
    tracknames.mask <- colnames(binstates)[mask]
    if (nrow(binstates) > 1) {
        statenames.sep <- apply(as.matrix(binstates[,mask]), 1, function(x) { tracknames.mask[x] })
        if (length(statenames.sep)==0) {
            stop("Something went wrong in constructing state names.")
        } else if (class(statenames.sep)=='list') {
            statenames <- sapply(statenames.sep, paste, collapse=sep)
        } else if (class(statenames.sep)=='matrix') {
            statenames <- apply(statenames.sep, 2, paste, collapse=sep)
        } else if (class(statenames.sep)=='character') {
            statenames <- statenames.sep
        }
    } else {
        statenames <- paste(tracknames.mask[binstates[,mask]], collapse=sep)
    }
    ## Convert to decimal
    decstates.all <- bin2dec(binstates)
    duplicate.mask <- !duplicated(decstates.all)
    decstates.all <- decstates.all[duplicate.mask]
    names(decstates.all) <- statenames[duplicate.mask]
    
    ### Select exclusive states
    if (!is.null(exclusive.table)) {
        if (is.null(tracks2compare) | is.null(conditions)) {
            stop("Arguments 'tracks2compare' and 'conditions' must be specified if 'exclusive.table' was specified.")
        }

        if (is.character(exclusive.table)) {
            excl.table <- utils::read.table(exclusive.table, header=TRUE, comment.char='#')
            if (!all(colnames(exclusive.table) == c('mark','group'))) {
                stop("Your 'exclusive.table' must be a tab-separated file with column names 'mark' and 'group'.")
            }
        } else if (is.data.frame(exclusive.table)) {
            excl.table <- exclusive.table
        } else {
            stop("Argument 'exclusive.table' must be a data.frame or a tab-separated file.")
        }
        
        loci <- split(excl.table$mark, excl.table$group)

        tracknames.split <- split(tracknames, conditions)
        # tracks2compare.split <- split(tracks2compare, conditions)
        for (cond in unique(conditions)) {
            names <- tracknames.split[[as.character(cond)]]
            
            for (locus in loci) {
                # indexes (no replicates) of histone modifications at locus
                track.index <- which((grepl(paste0(unlist(locus), collapse='|'), names)) & !duplicated(names))

                if (length(track.index) > 1) { # no restrictions if only one modification exists
                    mask <- rowSums(as.matrix(binstates[,names[track.index]]), na.rm = TRUE) < 2
                    binstates <- binstates[mask,]
                    
                    if (class(binstates)!='matrix') {
                        binstates <- matrix(binstates, ncol=numtracks)
                        colnames(binstates) <- tracknames
                    }
                }
            }
        }
    }
    
    ### Select specified differential states ###
    if (!is.null(diffstatespec)) {
        diffgroups <- levels(factor(diffstatespec))
        for (diffgroup in diffgroups) {
            track.index <- which(diffstatespec==diffgroup)
            if (grepl('^d\\.', diffgroup)) {
                mask0 <- !apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('|', x) })
                mask1 <- apply(as.matrix(binstates[,track.index]), 1, function(x) { Reduce('&', x) })
                mask <- !(mask0 | mask1)
            } else if (grepl('^x\\.', diffgroup)) {
                mask <- rep(TRUE, nrow(binstates))
            }
            binstates <- binstates[mask,]
            if (class(binstates)!='matrix') {
                binstates <- matrix(binstates, ncol=length(binstates))
                colnames(binstates) <- tracknames
            }
        }
    }

    ### Select differential states by min.diff ###
    # Idea: Get previously computed binary states and select those that are differential.
    # Beware: The column order is different from above
    binstates.diff <- NULL
    if (differential.states) {
        if (is.null(tracks2compare)) {
            stop("argument 'tracks2compare' must be specified if 'conditions' was specified")
        }
        tracknames.conditions <- unlist(split(tracknames, conditions))
        tracks2compare.split <- split(tracks2compare, conditions)
        intersect.tracks <- Reduce(intersect, lapply(tracks2compare.split, unique))
        bindiffmatrix <- dec2bin(0:(2^length(intersect.tracks)-1)) # all possible binary combinations (rows) between tracks (cols)
        controlsum <- apply(bindiffmatrix, 1, sum)
        bindiffmatrix <- bindiffmatrix[controlsum >= min.diff,]
        if (class(bindiffmatrix)!='matrix') { # R-behaviour differs with only one column
            bindiffmatrix <- matrix(bindiffmatrix, nrow=1)
        }
        colnames(bindiffmatrix) <- intersect.tracks
        diffstatespec.list <- list()
        ## Go through conditions
        for (condition in unique(conditions)) {
            tracks <- tracks2compare.split[[as.character(condition)]]
            #TODO: tracksNOT2use
            tracks2use <- tracks[tracks %in% intersect.tracks]
            num.tracks.split <- length(tracks2use)
            ## Replicate-expanded (col) bindiffmatrix
            bindiffmatrix.reps <- matrix(NA, nrow=nrow(bindiffmatrix), ncol=length(tracks2use), dimnames=list(combination=1:nrow(bindiffmatrix), tracks2use=tracks2use))
            for (track in intersect.tracks) { # go through replicates
                index <- which(track==tracks2use)
                bindiffmatrix.reps[,index] <- rep(bindiffmatrix[,as.character(track)], length(index))
            }
            ## Replace FALSE and TRUE by "0." and "1."
            if (ncol(bindiffmatrix.reps)==1) { # R-behaviour differs with only one column
                diffstatespec.reps <- t(t(apply(bindiffmatrix.reps, 1, function(x) { c('x.','d.')[x+1] })))
            } else {
                diffstatespec.reps <- t(apply(bindiffmatrix.reps, 1, function(x) { c('x.','d.')[x+1] }))
            }
            colnames(diffstatespec.reps) <- colnames(bindiffmatrix.reps)
            ## Attach the trackname to "0." or "1." for the next step
            if (ncol(diffstatespec.reps)==1) { # R-behaviour differs with only one column
                diffstatespec.list[[as.character(condition)]] <- t(t(apply(diffstatespec.reps, 1, function(x) { paste0(x, tracks2use) })))
            } else {
                diffstatespec.list[[as.character(condition)]] <- t(apply(diffstatespec.reps, 1, function(x) { paste0(x, tracks2use) }))
            }
            colnames(diffstatespec.list[[as.character(condition)]]) <- tracks2use
        }
        diffstatespecs <- do.call(cbind, diffstatespec.list)

        ## Go through all diffstate specifications
        binstates.list <- list()
        for (irow in 1:nrow(diffstatespecs)) {
            diffstatespec <- diffstatespecs[irow,]
            diffgroups <- levels(factor(diffstatespec))
            binstates.irow <- binstates[,tracknames.conditions] # reorder to match currently used ordering scheme
            for (diffgroup in diffgroups) {
                track.index <- which(diffstatespec==diffgroup)
                if (grepl('^d\\.', diffgroup)) {
                    mask0 <- !apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('|', x) }) # rows where all group members are 0
                    mask1 <- apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('&', x) }) # rows where all group members are 1
                    mask <- !(mask0 | mask1) # rows where not all group members are either 0 or 1
                } else if (grepl('^x\\.', diffgroup)) {
                    mask <- rep(TRUE, nrow(binstates.irow))
                }
                binstates.irow <- binstates.irow[mask,]
                if (class(binstates.irow)!='matrix') {
                    binstates.irow <- matrix(binstates.irow, ncol=length(binstates.irow))
                    colnames(binstates.irow) <- tracknames.conditions
                }
            }
            binstates.list[[irow]] <- binstates.irow
        }
        binstates.diff <- do.call(rbind, binstates.list)
        # Reorder to match initial order
        binstates.diff <- binstates.diff[,tracknames]
    }
    # There are still duplicate rows at this point, they are removed at the end

    ### Select common states ###
    # Idea: Create common states from scratch.
    # Beware: The column order is different from the binary states in the beginning but the same as in differential state selection by min.diff
    binstates.common <- NULL
    if (common.states) {
        if (is.null(tracks2compare)) {
            stop("argument 'tracks2compare' must be specified if 'conditions' was specified")
        }
        tracknames.conditions <- unlist(split(tracknames, conditions))
        tracks2compare.split <- split(tracks2compare, conditions)
        intersect.tracks <- Reduce(intersect, lapply(tracks2compare.split, unique))
        bincommonmatrix <- dec2bin(0:(2^length(intersect.tracks)-1)) # all possible binary combinations (rows) between tracks (cols)
        if (class(bincommonmatrix)!='matrix') { # R-behaviour differs with only one column
            bincommonmatrix <- matrix(bincommonmatrix, nrow=1)
        }
        colnames(bincommonmatrix) <- intersect.tracks
        commonstatespec.list <- list()
        ## Go through conditions
        for (condition in unique(conditions)) {
            tracks <- tracks2compare.split[[as.character(condition)]]
            #TODO: tracksNOT2use
            tracks2use <- tracks[tracks %in% intersect.tracks]
            num.tracks.split <- length(tracks2use)
            ## Replicate-expanded (col) bincommonmatrix
            bincommonmatrix.reps <- matrix(NA, nrow=nrow(bincommonmatrix), ncol=length(tracks2use), dimnames=list(combination=1:nrow(bincommonmatrix), tracks2use=tracks2use))
            for (track in intersect.tracks) { # go through replicates
                index <- which(track==tracks2use)
                bincommonmatrix.reps[,index] <- rep(bincommonmatrix[,as.character(track)], length(index))
            }
            ## Replace FALSE and TRUE by "0." and "1."
            if (ncol(bincommonmatrix.reps)==1) { # R-behaviour differs with only one column
                commonstatespec.reps <- t(t(apply(bincommonmatrix.reps, 1, function(x) { c('0.','1.')[x+1] })))
            } else {
                commonstatespec.reps <- t(apply(bincommonmatrix.reps, 1, function(x) { c('0.','1.')[x+1] }))
            }
            colnames(commonstatespec.reps) <- colnames(bincommonmatrix.reps)
            ## Attach the trackname to "0." or "1." for the next step
            if (ncol(commonstatespec.reps)==1) { # R-behaviour differs with only one column
                commonstatespec.list[[as.character(condition)]] <- t(t(apply(commonstatespec.reps, 1, function(x) { paste0(x, tracks2use) })))
            } else {
                commonstatespec.list[[as.character(condition)]] <- t(apply(commonstatespec.reps, 1, function(x) { paste0(x, tracks2use) }))
            }
            colnames(commonstatespec.list[[as.character(condition)]]) <- tracks2use
        }
        commonstatespecs <- do.call(cbind, commonstatespec.list)

        ## Go through all commonstate specifications
        binstates.list <- list()
        for (irow in 1:nrow(commonstatespecs)) {
            commonstatespec <- commonstatespecs[irow,]
            commongroups <- levels(factor(commonstatespec))
            binstates.irow <- binstates
            for (commongroup in commongroups) {
                track.index <- which(commonstatespec==commongroup)
                if (grepl('^0\\.', commongroup)) {
                    mask <- !apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('|', x) })
                } else if (grepl('^1\\.', commongroup)) {
                    mask <- apply(as.matrix(binstates.irow[,track.index]), 1, function(x) { Reduce('&', x) })
                }
                binstates.irow <- binstates.irow[mask,]
                if (class(binstates.irow)!='matrix') {
                    binstates.irow <- matrix(binstates.irow, ncol=length(binstates.irow))
                    colnames(binstates.irow) <- tracknames.conditions
                }
            }
            binstates.list[[irow]] <- binstates.irow
        }
        binstates.common <- do.call(rbind, binstates.list)
        # Reorder to match initial order
        binstates.common <- binstates.common[,tracknames]
    }
        
    ## Merge common and differential states
    if (differential.states | common.states) {
        binstates <- rbind(binstates.diff, binstates.common)
    }

    ### Construct state names ###
    mask <- !grepl('^r\\.', statespec) | !duplicated(statespec) # only first replicate
    tracknames.mask <- colnames(binstates)[mask]
    if (nrow(binstates) > 1) {
        statenames.sep <- apply(as.matrix(binstates[,mask]), 1, function(x) { tracknames.mask[x] })
        if (length(statenames.sep)==0) {
            stop("Something went wrong in constructing state names.")
        } else if (class(statenames.sep)=='list') {
            statenames <- sapply(statenames.sep, paste, collapse=sep)
        } else if (class(statenames.sep)=='matrix') {
            statenames <- apply(statenames.sep, 2, paste, collapse=sep)
        } else if (class(statenames.sep)=='character') {
            statenames <- statenames.sep
        }
    } else {
        statenames <- paste(tracknames.mask[binstates[,mask]], collapse=sep)
    }
    ## Convert to decimal
    decstates <- bin2dec(binstates)
    duplicate.mask <- !duplicated(decstates)
    decstates <- decstates[duplicate.mask]
    names(decstates) <- statenames[duplicate.mask]

    ## Put in data.frame and assign levels
    if (length(decstates)>0) {
        dec.df <- data.frame(combination=factor(paste0('[',names(decstates),']'), levels=paste0('[',names(decstates.all),']')), state=factor(decstates, levels=decstates.all))
        rownames(dec.df) <- NULL
    } else {
        dec.df <- data.frame(combination=factor(), state=integer())
    }

    return(dec.df)

}
