

#' Read bed-file into GRanges
#'
#' This is a simple convenience function to read a bed(.gz)-file into a \code{\link{GRanges}} object. The bed-file is expected to have the following fields: \code{chromosome, start, end, name, score, strand}.
#'
#' @param bedfile Filename of the bed or bed.gz file.
#' @param col.names A character vector giving the names of the columns in the \code{bedfile}. Must contain at least \code{c('chromosome','start','end')}.
#' @param col.classes A character vector giving the classes of the columns in \code{bedfile}. Speeds up the import.
#' @param skip Number of lines to skip at the beginning.
#' @param chromosome.format Desired format of the chromosomes. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).
#' @return A \code{\link{GRanges}} object with the contents of the bed-file.
#' @author Aaron Taudt
#' @importFrom utils read.table
#' @export
#'
#'@examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "liver-H3K4me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Import the file and skip the first 10 lines
#'data <- readCustomBedFile(bedfile, skip=10)
#'
readCustomBedFile <- function(bedfile, col.names=c('chromosome','start','end','name','score','strand'), col.classes=NULL, skip=0, chromosome.format='NCBI') {

    if (!is.null(col.classes)) {
        names(col.classes) <- col.names
        data <- utils::read.table(bedfile, col.names=col.names, colClasses=col.classes, skip=skip)
    } else {
        data <- utils::read.table(bedfile, col.names=col.names, skip=skip)
    }
    if (!is.null(data$strand)) {
        # GRanges compatible strand information
        data$strand <- sub('.','*',data$strand)
    } else {
        data$strand <- '*'
    }
    if (nrow(data) > 0) {
        # Adjust chromosome format
        data$chromosome <- sub('^chr', '', data$chromosome)
        if (chromosome.format=='UCSC') {
            data$chromosome <- paste0('chr', data$chromosome)
        }
    }
    # Convert to GRanges object
    gr <- methods::as(data, 'GRanges')
    start(gr) <- start(gr) + 1

    return(gr)

}

