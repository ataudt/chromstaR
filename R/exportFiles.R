#' Export genome browser uploadable files
#' 
#' These functions allow to export \code{\link{chromstaR-objects}} as files which can be uploaded to a genome browser. Peak calls are exported in BED format (.bed.gz), read counts in wiggle format (.wig.gz) as RPKM values, and combinatorial states are exported in BED format (.bed.gz).
#' 
#' @name exportFiles
#' @param model A \code{\link{chromstaR-objects}}.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either "_peaks.bed.gz" for peak-calls or "_counts.wig.gz" for read counts or "_combinations.bed.gz" for combinatorial states. Any existing file will be overwritten.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param separate.files A logical indicating whether or not to produce separate files for each track.
#' @param trackname Name that will be used in the "track name" field of the BED file.
#' @return \code{NULL}
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","combined_mode-differential.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'## Export peak calls and combinatorial states
#'exportPeaks(model, filename=tempfile())
#'exportCombinations(model, filename=tempfile())
#'
NULL


#' @describeIn exportFiles Export peak calls in BED format.
#' @export
exportPeaks <- function(model, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {
  
    model <- loadHmmsFromFiles(model)[[1]]
    if (class(model) == class.univariate.hmm) {
        exportUnivariatePeaks(list(model), filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if (class(model) == class.multivariate.hmm) {
        exportMultivariatePeaks(model, filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if (class(model) == class.combined.multivariate.hmm) {
        exportCombinedMultivariatePeaks(model, filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, trackname=trackname)
    }
  
}


#' @describeIn exportFiles Export read counts as RPKM values in wiggle format.
#' @export
exportCounts <- function(model, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {
  
    model <- loadHmmsFromFiles(model)[[1]]
    if (class(model) == 'GRanges') {
        exportBinnedData(list(model), filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if (class(model) == class.univariate.hmm) {
        exportUnivariateCounts(list(model), filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if (class(model) == class.multivariate.hmm) {
        exportMultivariateCounts(model, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if (class(model) == class.combined.multivariate.hmm) {
        exportCombinedMultivariateCounts(model, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
  
}


#' @describeIn exportFiles Export combinatorial states in BED format.
#' @param exclude.states A character vector with combinatorial states that will be excluded from export.
#' @param include.states A character vector with combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
#' @export
exportCombinations <- function(model, filename, header=TRUE, separate.files=TRUE, trackname=NULL, exclude.states='[]', include.states=NULL) {
  
    model <- loadHmmsFromFiles(model)[[1]]
    if (class(model) == class.multivariate.hmm) {
        exportMultivariateCombinations(model, filename=paste0(filename, '_combinations'), header=header, trackname=trackname, exclude.states=exclude.states, include.states=include.states)
    }
    if (class(model) == class.combined.multivariate.hmm) {
        exportCombinedMultivariateCombinations(model, filename=paste0(filename, '_combinations'), header=header, separate.files=separate.files, trackname=trackname, exclude.states=exclude.states, include.states=include.states)
    }
  
}