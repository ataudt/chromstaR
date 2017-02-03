#' Print uniHMM object
#' 
#' @param x An \code{\link{uniHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.uniHMM <- function(x, ...) {
    
    message("$info")
    print(x$info)
    message("\n$peaks")
    print(x$peaks)
    message("\nUse the list operator $ to access all elements of this object.")
  
}


#' Print multiHMM object
#' 
#' @param x An \code{\link{multiHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.multiHMM <- function(x, ...) {
    
    message("$info")
    print(x$info)
    message("\n$segments")
    print(x$segments)
    message("\nUse the list operator $ to access all elements of this object.")
  
}


#' Print combinedMultiHMM object
#' 
#' @param x An \code{\link{combinedMultiHMM}} object.
#' @param ... Ignored.
#' @return An invisible \code{NULL}.
#' @export
print.combinedMultiHMM <- function(x, ...) {
    
    message("$info")
    print(x$info)
    message("\n$segments")
    print(x$segments)
    message("\nUse the list operator $ to access all elements of this object.")
  
}
