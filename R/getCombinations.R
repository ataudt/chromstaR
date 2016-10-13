#' Get combinations
#' 
#' Get a DataFrame with combinations from a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param gr A \code{\link[GenomicRanges]{GRanges}} object from which the meta-data columns containing combinations will be extracted.
#' @return A DataFrame.
#' 
#' @export
#' @examples 
#'### Get an example multiHMM ###
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'### Get the combinations
#'bin.combs <- getCombinations(model$bins)
#'print(bin.combs)
#'seg.combs <- getCombinations(model$segments)
#'print(seg.combs)
#' 
getCombinations <- function(gr) {
    
    combs <- mcols(gr)[grepl('combination', names(mcols(gr)))]
    return(combs)
    
}