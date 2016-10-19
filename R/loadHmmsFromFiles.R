#' Load \pkg{chromstaR} objects from file
#'
#' Wrapper to load \pkg{\link{chromstaR}} objects from file and check the class of the loaded objects.
#'
#' @param files A list of \code{\link{chromstaR-objects}} or a vector of files that contain such objects.
#' @param check.class Any combination of \code{c('GRanges', 'uniHMM', 'multiHMM', 'combinedMultiHMM')}. If any of the loaded objects does not belong to the specified class, an error is thrown.
#' @return A list of \code{\link[chromstaR:chromstaR-objects]{chromstaR-object}}.
#' @export
#' @examples
#'## Get an example BAM file
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
#'temp.file <- tempfile()
#'save(hmm, file=temp.file)
#'loaded.hmm <- loadHmmsFromFiles(temp.file)[[1]]
#'class(loaded.hmm)
#'
loadHmmsFromFiles <- function(files, check.class=c('GRanges', 'uniHMM', 'multiHMM', 'combinedMultiHMM')) {

#     ptm <- startTimedMessage("Loading data from files ...")
    if (any(! check.class %in% c('GRanges', class.univariate.hmm, class.multivariate.hmm, class.combined.multivariate.hmm))) {
        stop("Argument 'check.class' must contain any combination of c('", paste0(c('GRanges', class.univariate.hmm, class.multivariate.hmm, class.combined.multivariate.hmm), collapse="', '"), "').")
    }
    modellist <- list()
    if (is.character(files)) {
        for (file in files) {
            temp.env <- new.env()
            model <- get(load(file, envir=temp.env), envir=temp.env)
            if (! class(model) %in% check.class) {
                stop("File '", file, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[basename(file)]] <- model
        }
    } else if (class(files) %in% check.class) {
        modellist[[1]] <- files
    } else if (is.list(files)) {
        for (file in files) {
            model <- file
            if (! class(model) %in% check.class) {
                stop("List entry '", length(modellist)+1, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[length(modellist)+1]] <- model
        }
        names(modellist) <- names(files)
    } else if (! class(files) %in% check.class) {
        stop("Input does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
    }
#     stopTimedMessage(ptm)
    return(modellist)

}
