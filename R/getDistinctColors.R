

#' Get distinct colors
#' 
#' Get a set of distinct colors selected from \code{\link{colors}}.
#' 
#' The function computes the euclidian distance between all \code{\link{colors}} and iteratively selects those that have the furthest closes distance to the set of already selected colors.
#' 
#' @param n Number of colors to select.
#' @param start.color Color to start the selection process from.
#' @param exclude.colors Character vector with colors that should not be used.
#' @param exclude.rgb.above Exclude colors where all RGB values are above. This is useful to exclude whitish colors.
#' @return A character vector with colors.
#' @author Aaron Taudt
#' @importFrom grDevices col2rgb
#' @importFrom stats dist
#' @export
#' @examples
#'cols <- getDistinctColors(5)
#'pie(rep(1,5), labels=cols, col=cols)
#'
getDistinctColors <- function(n, start.color='blue4', exclude.colors=c('white','black','gray','grey'), exclude.rgb.above=210) {
	
	cols <- grDevices::colors()
	# Exclude unwanted colors
	cols <- grep(paste(exclude.colors, collapse='|'), cols, invert=TRUE, value=TRUE)
	# Get RGB values
	rgbs <- t(grDevices::col2rgb(cols))
	rownames(rgbs) <- cols
	# Exclude whitish colors
	rgbs <- rgbs[apply(rgbs, 1, function(x) { !all(x>exclude.rgb.above) }), ]
	# Calculate distance
	coldist <- as.matrix(stats::dist(rgbs, method='euclidean'))
	
	# Iteratively select colors
	selected.cols <- character()
	selected.cols[1] <- start.color
	for (i1 in 2:n) {
		m <- as.matrix(coldist[,selected.cols])
		closest.dist <- apply(m, 1, min)
		furthest.dist <- which.max(closest.dist)
		selected.cols[i1] <- names(furthest.dist)
	}
	return(selected.cols)
	
}
