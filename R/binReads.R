

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed(.gz) format into read counts in equidistant windows.
#'
#' Convert aligned reads from .bam or .bed(.gz) files into read counts in equidistant windows (bins). This function uses \code{\link[GenomicRanges]{countOverlaps}} to calculate the read counts, or alternatively \code{\link[bamsignals]{bamProfile}} if option \code{use.bamsignals} is set (only effective for .bam files).
#'
#' @aliases binning
#' @param file A file with aligned reads. Alternatively a \code{\link{GRanges}} with aligned reads.
#' @param experiment.table An \code{\link{experiment.table}} containing the supplied \code{file}. This is necessary to uniquely identify the file in later steps of the workflow. Set to \code{NULL} if you don't have it (not recommended).
#' @inheritParams readBamFileAsGRanges
#' @inheritParams readBedFileAsGRanges
#' @param binsizes An integer vector specifying the bin sizes to use.
#' @param stepsizes An integer vector specifying the step size. One number can be given for each element in \code{binsizes}, \code{reads.per.bin} and \code{bins} (in that order).
#' @param bins A \code{\link{GRanges}} or a named \code{list()} with \code{\link{GRanges}} containing precalculated bins produced by \code{\link{fixedWidthBins}} or \code{\link{variableWidthBins}}. Names of the list must correspond to the binsize. If the list is unnamed, an attempt is made to automatically determine the binsize.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin size will be selected accordingly.
#' @param variable.width.reference A BAM file that is used as reference to produce variable width bins. See \code{\link{variableWidthBins}} for details.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param use.bamsignals If \code{TRUE} the \pkg{\link[bamsignals]{bamsignals}} package is used for parsing of BAM files. This gives tremendous speed advantage for only one binsize but linearly increases for multiple binsizes, while \code{use.bamsignals=FALSE} has a binsize dependent runtime and might be faster if many binsizes are calculated.
#' @param format One of \code{c('bed','bam','GRanges')} or \code{NULL} if the format is to be determined automatically.
#' @return If only one bin size was specified for option \code{binsizes}, the function returns a single \code{\link{GRanges}} object with meta data column 'counts' that contains the read count. If multiple \code{binsizes} were specified , the function returns a named \code{list()} of \link{GRanges} objects.
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom bamsignals bamCount
#' @export
#'
#'@examples
#'## Get an example BAM file with ChIP-seq reads
#'file <- system.file("extdata", "euratrans",
#'                    "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                     package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'data(experiment_table)
#'binned <- binReads(file, experiment.table=experiment_table,
#'                   assembly=rn4_chrominfo, binsizes=1000,
#'                   stepsizes=500, chromosomes='chr12')
#'print(binned)
#'
binReads <- function(file, experiment.table=NULL, assembly, bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000, blacklist=NULL, binsizes=1000, stepsizes=1/5 * binsizes, reads.per.bin=NULL, bins=NULL, variable.width.reference=NULL, use.bamsignals=TRUE, format=NULL) {

    ## Determine format
    if (is.null(format)) {
        if (is.character(file)) {
            file.clean <- sub('\\.gz$','', file)
            format <- rev(strsplit(file.clean, '\\.')[[1]])[1]
        } else if (class(file)=='GRanges') {
            format <- 'GRanges'
        } else {
            stop("Could not determine format automatically. Please specify it via the 'format' parameter.")
        }
    }
    if (! format %in% c('bed','bam','GRanges')) {
        stop("Unknown format. Needs to be one of 'bed', 'bam' or 'GRanges'.")
    }

    ## Sanity checks
    if (format=='bed') {
        temp <- assembly # trigger error if not defined
    }
    if (class(bins) == 'GRanges') {
        bins <- list(bins)
        names(bins) <- width(bins[[1]])[1]
    }
    if (class(bins) == 'list') {
        if (is.null(names(bins))) {
            names(bins) <- sapply(bins, function(x) { width(x)[1] })
        }
    }

    ## Variables for bamsignals
    paired.end <- 'ignore'
    if (pairedEndReads) {
        paired.end <- 'filter'
    }
    
    ## Create INFO object as row from the experiment.table
    if (!is.null(experiment.table)) {
        check.experiment.table(experiment.table)
        info <- experiment.table[basename(as.character(experiment.table$file))==basename(file),]
        ID <- paste0(info$mark, '-', info$condition, '-rep', info$rep)
        info$ID <- ID
        if (pairedEndReads != info$pairedEndReads) {
            warning("Option 'pairedEndReads' overwritten by entry in 'experiment.table'.")
        }
    } else {
        info <- NULL
        if (format=='GRanges') {
            ID <- 'GRanges'
        } else {
            ID <- basename(file)
        }
    }

    ### Read in the data
    data <- NULL
    if (format == "bed") {
        ## BED (0-based)
        if (!remove.duplicate.reads) {
            data <- readBedFileAsGRanges(file, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=FALSE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
        } else {
            data <- readBedFileAsGRanges(file, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=TRUE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
        }
        chrom.lengths <- seqlengths(data)
    } else if (format == "bam") {
        ## BAM (1-based)
        if (use.bamsignals) {
            ## Check if bamindex exists
            bamindex.raw <- sub('\\.bai$', '', bamindex)
            bamindex <- paste0(bamindex.raw,'.bai')
            if (!file.exists(bamindex)) {
                ptm <- startTimedMessage("Making bam-index file ...")
                bamindex.own <- Rsamtools::indexBam(file)
                # warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
                bamindex <- bamindex.own
                stopTimedMessage(ptm)
            }
            ptm <- startTimedMessage(paste0("Reading header from ", file, " ..."))
            chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(file))
            stopTimedMessage(ptm)
        } else {
            if (!remove.duplicate.reads) {
                data <- readBamFileAsGRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=FALSE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
            } else {
                data <- readBamFileAsGRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=TRUE, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
            }
            chrom.lengths <- seqlengths(data)
        }
    } else if (format == "GRanges") {
        ## GRanges (1-based)
        data <- file
        chrom.lengths <- seqlengths(data)
    } else {
        stop("Unknown file format.")
    }

    ## Select chromosomes to bin
    if (is.null(chromosomes)) {
        chromosomes <- names(chrom.lengths)
    }
    chroms2use <- intersect(chromosomes, names(chrom.lengths))
    ## Stop if none of the specified chromosomes exist
    if (length(chroms2use)==0) {
      chrstring <- paste0(chromosomes, collapse=', ')
      stop('Could not find length information for any of the specified chromosomes: ', chrstring)
    }
    skipped.chroms <- setdiff(chromosomes, chroms2use)
    if (length(skipped.chroms)>0) {
        warning("Could not find chromosomes ", paste0(skipped.chroms, collapse=', '), ".")
    }

    ## Select only desired chromosomes
    if (!is.null(data)) {
        ptm <- startTimedMessage("Subsetting specified chromosomes ...")
        data <- data[seqnames(data) %in% chroms2use]
        data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
        ## Drop seqlevels where seqlength is NA
        na.seqlevels <- seqlevels(data)[is.na(seqlengths(data))]
        data <- data[seqnames(data) %in% seqlevels(data)[!is.na(seqlengths(data))]]
        data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
        if (length(na.seqlevels) > 0) {
            warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
        }
        stopTimedMessage(ptm)
    }
  
    ### Make variable width bins ###
    if (!is.null(variable.width.reference)) {
        message("Making variable width bins:")
        if (is.character(variable.width.reference)) {
            variable.width.reference.clean <- sub('\\.gz$','', variable.width.reference)
            vformat <- rev(strsplit(variable.width.reference.clean, '\\.')[[1]])[1]
        } else if (class(variable.width.reference)=='GRanges') {
            vformat <- 'GRanges'
        }
        if (vformat == 'bam') {
            refreads <- readBamFileAsGRanges(variable.width.reference, bamindex=variable.width.reference, chromosomes=chroms2use, pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
        } else if (vformat == 'bed') {
            refreads <- readBedFileAsGRanges(variable.width.reference, assembly=assembly, chromosomes=chroms2use, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width, blacklist=blacklist)
        }
        bins.binsize <- variableWidthBins(refreads, binsizes=binsizes, chromosomes=chroms2use)
        message("Finished making variable width bins.")
    }
        
    ### Make fixed width bins ###
    if (is.null(variable.width.reference)) {
        bins.binsize <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes)
    }
  
    ### Make reads.per.bin bins ###
    bins.rpb <- NULL
    if (!is.null(data)) {
        numcountsperbp <- length(data) / sum(as.numeric(seqlengths(data)))
        binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
        bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb)
    } else if (use.bamsignals & !is.null(reads.per.bin)) {
        ptm <- startTimedMessage("Parsing bamfile to determine binsize for reads.per.bin option ...")
        bins.helper <- suppressMessages( fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=1e6)[[1]] )
        counts <- tryCatch({
            counts <- bamsignals::bamCount(file, bins.helper, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE)
        }, error = function(err) {
            counts <- bamsignals::bamCount(file, bins.helper, mapqual=min.mapq, paired.end=paired.end, paired.end.max.frag.length=max.fragment.width, verbose=FALSE)
        })
        stopTimedMessage(ptm)
        numcountsperbp <- sum(as.numeric(counts)) / sum(as.numeric(chrom.lengths[chroms2use]))
        binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
        bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb)
    }
    
    ### Combine in bins.list ###
    bins.list <- c(bins.binsize, bins.rpb, bins)
    if (length(stepsizes) != length(bins.list)) {
        stop("Need one step size for each binsize. Make sure that length(stepsizes) == length(binsizes) + length(reads.per.bin) + length(bins).")
    }
 
    ### Loop over all binsizes ###
    if (!use.bamsignals | format=='bed' | format=='GRanges') {
        ptm <- startTimedMessage("Splitting into strands ...")
        data.plus <- data[strand(data)=='+']
        data.minus <- data[strand(data)=='-']
        data.star <- data[strand(data)=='*']
        ptm <- stopTimedMessage(ptm)
    }
    for (ibinsize in 1:length(bins.list)) {
        binsize <- as.numeric(names(bins.list)[ibinsize])
        bins <- bins.list[[ibinsize]]
        stepsize <- stepsizes[ibinsize]
        numsteps <- ceiling(binsize / stepsize)
        offsets <- stepsize * ((1:numsteps) - 1)
        acounts <- array(NA, dim = c(length(bins), numsteps), dimnames = list(bin=NULL, offset=offsets))
        for (ioffset in 1:numsteps) {
            offset <- offsets[ioffset]
            bins.offset <- suppressWarnings( shift(bins, shift = offset) )
        
            if (format == 'bam' & use.bamsignals) {
                ptm <- startTimedMessage("Counting overlaps for binsize ", binsize, " with offset ", offset, " ...")
                bins.offset$counts <- tryCatch({
                    bins.offset$counts <- bamsignals::bamCount(file, bins.offset, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE)
                }, error = function(err) {
                    bins.offset$counts <<- bamsignals::bamCount(file, bins.offset, mapqual=min.mapq, paired.end=paired.end, paired.end.max.frag.length=max.fragment.width, verbose=FALSE)
                })
                stopTimedMessage(ptm)
                
            } else {
                readsperbin <- round(length(data) / sum(as.numeric(seqlengths(data))) * binsize, 2)
                ptm <- startTimedMessage("Counting overlaps for binsize ",binsize," with offset ", offset, ", with on average ",readsperbin," reads per bin ...")
                scounts <- suppressWarnings( GenomicRanges::countOverlaps(bins.offset, data.star) )
                mcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins.offset, data.minus) )
                pcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins.offset, data.plus) )
                counts <- mcounts + pcounts + scounts
                countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
                colnames(countmatrix) <- c('counts','mcounts','pcounts')
                bins.offset$counts <- countmatrix[, 'counts']
                stopTimedMessage(ptm)
        
            }
            acounts[, as.character(offset)] <- mcols(bins.offset)[,'counts']
        
        }
        
        bins <- bins.list[[ibinsize]]
        bins$counts <- acounts
        bins.list[[ibinsize]] <- bins
        attr(bins.list[[ibinsize]], 'info') <- info
        attr(bins.list[[ibinsize]], 'min.mapq') <- min.mapq
            
    } ### end loop binsizes ###

    if (length(bins.list) == 1) {
        return(bins.list[[1]])
    } else {
        return(bins.list)
    }


}


