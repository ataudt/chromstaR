
#=====================================================
# Helper functions   
#=====================================================
insertchr <- function(gr) {
    # Change chromosome names from '1' to 'chr1' if necessary
    mask <- which(!grepl('chr', seqnames(gr)))
    mcols(gr)$chromosome <- as.character(seqnames(gr))
    mcols(gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(gr)$chromosome[mask])
    mcols(gr)$chromosome <- as.factor(mcols(gr)$chromosome)
    return(gr)
}

#=====================================================
# Export binned data
#=====================================================
# #' Export genome browser viewable files
# #'
# #' Export read counts (RPKM) as genome browser viewable file
# #'
# #' Export read counts from \code{\link{binned.data}} as a file which can be uploaded into a genome browser. Read counts are exported in WIGGLE format as RPKM (.wig.gz).
# #'
# #' @author Aaron Taudt
# #' @param binned.data.list A \code{list()} of \code{\link{binned.data}} objects or vector of files that contain such objects.
# #' @param filename The name of the file that will be written. The ending ".wig.gz" for read counts will be appended. Any existing file will be overwritten.
# #' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
# #' @param separate.files A logical indicating whether or not to produce separate files for each object in \code{binned.data.list}.
# #' @param trackname Name that will be used in the "track name" field of the file.
# #' @return \code{NULL}
# #' @seealso \code{\link{exportUnivariates}}, \code{\link{exportMultivariate}}
# #' @importFrom utils write.table
# #' @importFrom grDevices col2rgb
# #' @examples
# #'## Get an example BAM file
# #'file <- system.file("extdata", "euratrans",
# #'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
# #'                        package="chromstaRData")
# #'## Bin the file into bin size 1000bp
# #'data(rn4_chrominfo)
# #'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
# #'                   chromosomes='chr12')
# #'## Export the binned read counts
# #'exportBinnedData(list(binned), filename=tempfile())
# #'
exportBinnedData <- function(binned.data.list, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {

    ## Load data
    binned.data.list <- loadHmmsFromFiles(binned.data.list)

    ## Transform to GRanges
    binned.data.list <- lapply(binned.data.list, insertchr)

    # Variables
    nummod <- length(binned.data.list)
    filename <- paste0(filename,".wig.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }
    readcol <- paste(grDevices::col2rgb(getStateColors('counts')), collapse=',')

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        message('Writing track ',imod,' / ',nummod)
        b <- binned.data.list[[imod]]
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', imod, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            message('Writing to file ',filename.sep)
            cat("", file=filename.gz)
        }
        priority <- 50 + 4*imod
        binsize <- width(b[1])
        name <- names(binned.data.list)[imod]
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("read count for ", name)
            } else {
                trackname.string <- paste0("read count for ", name, ", ", trackname)
            }
            cat(paste0('track type=wiggle_0 name="',trackname.string,'" description="',trackname.string,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        b$counts <- rpkm.vector(b$counts, binsize=mean(width(b))) # RPKM normalization
        for (chrom in unique(b$chromosome)) {
            cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(b[b$chromosome==chrom]$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
        }
    }
    if (!separate.files) {
        close(filename.gz)
    }
}


#=====================================================
# Export univariate HMMs
#=====================================================
# #' Export genome browser viewable files
# #'
# #' Export univariate peak-calls and read counts (RPKM) as genome browser viewable file
# #'
# #' Export \code{\link{uniHMM}} objects as files which can be uploaded into a genome browser. Peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format as RPKM (.wig.gz).
# #'
# #' @author Aaron Taudt
# #' @param hmm.list A \code{list()} of \code{\link{uniHMM}} objects or vector of files that contain such objects.
# #' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for peak-calls or ".wig.gz" for read counts. Any existing file will be overwritten.
# #' @param what A character vector specifying what will be exported. Supported are \code{c('peaks', 'counts')}.
# #' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
# #' @param separate.files A logical indicating whether or not to produce separate files for each hmm in \code{hmm.list}.
# #' @param trackname Name that will be used in the "track name" field of the file.
# #' @return \code{NULL}
# #' @seealso \code{\link{exportBinnedData}}, \code{\link{exportMultivariate}}
# #' @examples
# #'## Get an example BAM file
# #'file <- system.file("extdata", "euratrans",
# #'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
# #'                        package="chromstaRData")
# #'## Bin the file into bin size 1000bp
# #'data(rn4_chrominfo)
# #'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
# #'                   chromosomes='chr12')
# #'## Fit the univariate Hidden Markov Model
# #'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
# #'## Export
# #'exportUnivariates(list(hmm), filename=tempfile(), what=c('peaks','counts'))
# #'
exportUnivariates <- function(hmm.list, filename, what=c('peaks', 'counts'), header=TRUE, separate.files=TRUE, trackname=NULL) {
    if ('peaks' %in% what) {
        exportUnivariatePeaks(hmm.list, filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if ('counts' %in% what) {
        exportUnivariateCounts(hmm.list, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
}

#----------------------------------------------------
# Export peak-calls from univariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportUnivariatePeaks <- function(hmm.list, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {

    ## Load models
    hmm.list <- loadHmmsFromFiles(hmm.list, check.class=class.univariate.hmm)

    ## Transform to GRanges
    peaklist <- lapply(hmm.list, '[[', 'peaks')
    peaklist <- lapply(peaklist, insertchr)

    # Variables
    nummod <- length(hmm.list)
    filename <- paste0(filename,".bed.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }

    # Write first line to file
    if (!separate.files) {
        ptm <- startTimedMessage('Writing to file ',filename, ' ...')
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        hmm <- hmm.list[[imod]]
        if (is.null(hmm$info$ID)) {
            ID <- paste0('track-',imod)
        } else {
            ID <- hmm$info$ID
        }
        if (separate.files) {
            filename.sep <- paste0(sub('.bed.gz$', '', filename), '_', ID, '.bed.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        }
        peaks <- peaklist[[imod]]
        priority <- 51 + 4*imod
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("univariate peak calls for ", ID)
            } else {
                trackname.string <- paste0("univariate peak calls for ", ID, ", ", trackname)
            }
            cat(paste0("track name=\"",trackname.string,"\" description=\"",trackname.string,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
        }
        if (is.null(peaks$peakScores)) {
            peaks$peakScores <- 0
        }
        df <- as.data.frame(peaks)
        df$peakNumber <- paste0('peak_', 1:nrow(df))
        df$strand <- sub('\\*', '.', df$strand)
        df <- df[,c('chromosome','start','end','peakNumber','peakScores','strand')]
        # Make score integer
        df$peakScores <- round(df$peakScores*1000)
        numsegments <- nrow(df)
        df <- cbind(df, thickStart=df$start, thickEnd=df$end)
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$thickStart - 1
        # Colors
        RGB <- t(grDevices::col2rgb(getStateColors('modified')))
        df$itemRgb <- apply(RGB,1,paste,collapse=",")
        if (nrow(df) == 0) {
            warning('hmm ',imod,' does not contain any \'modified\' calls')
        } else {
            utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
            stopTimedMessage(ptm)
        }
    }
    if (!separate.files) {
        close(filename.gz)
        stopTimedMessage(ptm)
    }

}


#----------------------------------------------------
# Export read counts from univariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportUnivariateCounts <- function(hmm.list, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {

    ## Load models
    hmm.list <- loadHmmsFromFiles(hmm.list, check.class=class.univariate.hmm)

    ## Transform to GRanges
    grl <- lapply(hmm.list, '[[', 'bins')
    hmm.grl <- lapply(grl, insertchr)

    # Variables
    nummod <- length(hmm.list)
    filename <- paste0(filename,".wig.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }
    readcol <- paste(grDevices::col2rgb(getStateColors('counts')), collapse=',')

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        hmm <- hmm.list[[imod]]
        if (is.null(hmm$info$ID)) {
            ID <- paste0('track-',imod)
        } else {
            ID <- hmm$info$ID
        }
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', ID, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('  Writing track ',imod,' / ',nummod, ' ...')
        }
        hmm.gr <- hmm.grl[[imod]]
        priority <- 50 + 4*imod
        binsize <- width(hmm.gr[1])
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("read count for ", ID)
            } else {
                trackname.string <- paste0("read count for ", ID, ", ", trackname)
            }
            cat(paste0('track type=wiggle_0 name="',trackname.string,'" description="',trackname.string,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        hmm.gr$counts <- rpkm.vector(hmm.gr$counts, binsize=mean(width(hmm.gr))) # RPKM normalization
        for (chrom in unique(hmm.gr$chromosome)) {
            cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(hmm.gr[hmm.gr$chromosome==chrom]$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
            stopTimedMessage(ptm)
        }
    }
    if (!separate.files) {
        close(filename.gz)
        stopTimedMessage(ptm)
    }
}


#=====================================================
# Export multivariate HMMs
#=====================================================
# #' Export genome browser viewable files
# #'
# #' Export multivariate calls and read counts (RPKM) as genome browser viewable file
# #'
# #' Export \code{\link{uniHMM}} objects as files which can be uploaded into a genome browser. Combinatorial states and peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format as RPKM (.wig.gz).
# #'
# #' @author Aaron Taudt
# #' @param hmm A \code{\link{multiHMM}} object or file that contains such an object.
# #' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for combinatorial states and peak-calls or ".wig.gz" for read counts. Any existing file will be overwritten.
# #' @param what A character vector specifying what will be exported. Supported are \code{c('combinations', 'peaks', 'counts')}.
# #' @param exclude.states A character vector with combinatorial states that will be excluded from export.
# #' @param include.states A character vector with combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
# #' @param trackname Name that will be used in the "track name" field of the BED file.
# #' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
# #' @param separate.files A logical indicating whether or not to produce separate files for peaks if \code{what} contains 'peaks' or 'counts'.
# #' @return \code{NULL}
# #' @seealso \code{\link{exportUnivariates}}, \code{\link{exportBinnedData}}
# #' @examples
# #'## Get an example multiHMM
# #'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
# #'                     package="chromstaR")
# #'model <- get(load(file))
# #'## Export peak calls and combinatorial states
# #'exportMultivariate(model, filename=tempfile(), what=c('peaks','combinations'))
# #'
exportMultivariate <- function(hmm, filename, what=c('combinations', 'peaks', 'counts'), exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=TRUE) {
    if ('combinations' %in% what) {
        exportMultivariateCombinations(hmm, filename=paste0(filename, '_combinations'), exclude.states=exclude.states, include.states=include.states, trackname=trackname, header=header)
    }
    if ('peaks' %in% what) {
        exportMultivariatePeaks(hmm, filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, trackname=trackname)
    }
    if ('counts' %in% what) {
        exportMultivariateCounts(hmm, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
}

#----------------------------------------------------
# Export combinatorial states or peak-calls from multivariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportMultivariateCombinations <- function(hmm, filename, separate.tracks=TRUE, exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE) {

    ## Intercept user input
    if (is.data.frame(exclude.states)) {
        exclude.states <- exclude.states$combination
    }
    if (is.data.frame(include.states)) {
        include.states <- include.states$combination
    }

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.multivariate.hmm)[[1]]

    ## Export the combinatorial states
    segments <- hmm$segments

    # Exclude and include states
    combstates <- levels(segments$combination)
    if (is.null(include.states)) {
        combstates2use <- setdiff(combstates, exclude.states)
    } else {
        combstates2use <- intersect(combstates, include.states)
    }
    segments <- segments[mcols(segments)$combination %in% combstates2use]
    segments <- insertchr(segments)
    segments.df <- as.data.frame(segments)

    # Make header
    if (is.null(trackname)) {
        trackname <- paste0('combinations')
    } else {
        trackname <- paste0('combinations, ', trackname)
    }

    ## Export combinatorial states
    # Generate the colors for each combinatorial state
    colors <- getDistinctColors(length(levels(segments.df$combination)))
    RGBs <- t(grDevices::col2rgb(colors))
    RGBlist <- list()
    for (i1 in 1:3) {
        RGBlist[[i1]] <- RGBs[,i1]
    }
    RGBlist$sep <- ','
    RGBs <- do.call(paste, RGBlist)
    itemRgb <- RGBs[as.integer(factor(segments.df$combination, levels=sort(levels(segments.df$combination))))]

    ## Write first line to file
    filename <- paste0(filename,".bed.gz")
    filename.gz <- gzfile(filename, 'w')
    ptm <- startTimedMessage('Writing to file ',filename, ' ...')
    cat("", file=filename.gz)

    # Write to file
    numsegments <- nrow(segments.df)
    if (is.null(segments.df$differential.score)) {
        segments.df$differential.score <- 0
    }
    names(segments.df)[names(segments.df)=='differential.score'] <- 'score'
    df <- cbind(segments.df[,c('chromosome','start','end','combination','score')], strand=rep(".",numsegments), thickStart=segments.df$start, thickEnd=segments.df$end, itemRgb=itemRgb)
    # Make score integer
    df$score <- round(df$score*1000)
    # Convert from 1-based closed to 0-based half open
    df$start <- df$start - 1
    df$thickStart <- df$thickStart - 1
    if (header) {
        cat(paste0('track name="',trackname, '" description="', trackname, '" visibility=1 itemRgb=On priority=49\n'), file=filename.gz, append=TRUE)
    }
    utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    
    close(filename.gz)
    stopTimedMessage(ptm)

}

#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportMultivariatePeaks <- function(hmm, filename, trackname=NULL, header=TRUE, separate.files=TRUE) {

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.multivariate.hmm)[[1]]

    ## Write first line to file
    if (!separate.files) {
        filename <- paste0(filename,".bed.gz")
        filename.gz <- gzfile(filename, 'w')
        cat("", file=filename.gz)
    }

    ## Export peaks on a per bin basis to obtain proper posterior scores
    for (imod in 1:length(hmm$info$ID)) {
        ID <- hmm$info$ID[imod]
        
        ## Select only segments with peaks
        peaks <- insertchr(hmm$peaks[[ID]])
        peaks.df <- as.data.frame(peaks)
        peaks.df$peakNumber <- paste0('peak_', 1:nrow(peaks.df))
        
        # Data.frame for write.table
        df <- peaks.df[,c('chromosome','start','end','peakNumber','peakScores','strand')]
        df$strand <- sub('\\*', '.', df$strand)
        
        # Make score integer
        df$peakScores <- round(df$peakScores*1000)
        
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$start
        df$thickEnd <- df$end
        # Colors
        RGB <- t(grDevices::col2rgb(getStateColors('modified')))
        df$itemRgb <- apply(RGB,1,paste,collapse=",")
        
        ## Write to file
        if (separate.files) {
            filename.sep <- paste0(sub('.bed.gz$', '', filename), '_', ID, '.bed.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('Writing to file ',filename, ' ...')
        }
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("peaks for ", ID)
            } else {
                trackname.string <- paste0("peaks for ", ID, ", ", trackname)
            }
            priority <- 52 + 4*imod
            cat(paste0('track name="', trackname.string, '" description="', trackname.string, '" visibility=1 itemRgb=On priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        if (separate.files) {
            close(filename.gz)
        }
        stopTimedMessage(ptm)
    }
    if (!separate.files) {
        close(filename.gz)
    }
}

#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportMultivariateCounts <- function(hmm, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.multivariate.hmm)[[1]]

    ## RPKM normalization
    hmm$bins$counts <- rpkm.matrix(hmm$bins$counts, binsize=mean(width(hmm$bins)))
    
    ## Variables
    filename <- paste0(filename,".wig.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }
    nummod <- length(hmm$info$ID)
    readcol <- paste(grDevices::col2rgb(getStateColors('counts')), collapse=',')

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        ID <- hmm$info$ID[imod]
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', ID, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('  Writing track ',imod,' / ',nummod, ' ...')
        }
        priority <- 50 + 4*imod
        binsize <- width(hmm$bins[1])
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("read count for ", ID)
            } else {
                trackname.string <- paste0("read count for ", ID, ", ", trackname)
            }
            cat(paste0('track type=wiggle_0 name="',trackname.string,'" description="',trackname.string,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        for (chrom in seqlevels(hmm$bins)) {
            if (!grepl('chr', chrom)) {
                chromr <- sub(pattern='^', replacement='chr', chrom)
            } else {
                chromr <- chrom
            }
            cat(paste0("fixedStep chrom=",chromr," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(hmm$bins[seqnames(hmm$bins)==chrom]$counts[,imod], file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
        }
        stopTimedMessage(ptm)
    }
    if (!separate.files) {
        close(filename.gz)
    }
}


#=====================================================
# Export combined multivariate HMMs
#=====================================================
# #' Export genome browser viewable files
# #'
# #' Export multivariate calls as genome browser viewable file
# #'
# #' Export \code{\link{combinedMultiHMM}} objects as files which can be uploaded into a genome browser. Combinatorial states are exported in BED format (.bed.gz). Read counts are exported in WIGGLE format as RPKM (.wig.gz).
# #'
# #' @author Aaron Taudt
# #' @param hmm A \code{\link{combinedMultiHMM}} object or file that contains such an object.
# #' @param filename The name of the file that will be written. The ending ".bed.gz" for combinatorial states will be appended. Any existing file will be overwritten.
# #' @param what A character vector specifying what will be exported. Supported are \code{c('combinations','peaks','counts')}.
# #' @param exclude.states A vector of combinatorial states that will be excluded from export.
# #' @param include.states A vector of combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
# #' @param trackname Name that will be used in the "track name" field of the BED file.
# #' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
# #' @param separate.files A logical indicating whether or not to produce separate files for each condition.
# #' @return \code{NULL}
# #' @examples
# #'## Get an example multiHMM
# #'file <- system.file("data","combined_mode-differential.RData",
# #'                     package="chromstaR")
# #'model <- get(load(file))
# #'## Export peak calls and combinatorial states
# #'exportCombinedMultivariate(model, filename=tempfile(), what=c('peaks','combinations'))
# #'
exportCombinedMultivariate <- function(hmm, filename, what=c('combinations','peaks'), exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=TRUE) {
    if ('combinations' %in% what) {
        exportCombinedMultivariateCombinations(hmm, filename=paste0(filename, '_combinations'), exclude.states=exclude.states, include.states=include.states, trackname=trackname, header=header, separate.files=separate.files)
    }
    if ('peaks' %in% what) {
        exportCombinedMultivariatePeaks(hmm, filename=paste0(filename, '_peaks'), trackname=trackname, header=header, separate.files=separate.files)
    }
    if ('counts' %in% what) {
        exportCombinedMultivariateCounts(hmm, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files, trackname=trackname)
    }
}

#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportCombinedMultivariatePeaks <- function(hmm, filename, trackname=NULL, header=TRUE, separate.files=TRUE) {

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.combined.multivariate.hmm)[[1]]

    ## Write first line to file
    if (!separate.files) {
        filename <- paste0(filename,".bed.gz")
        filename.gz <- gzfile(filename, 'w')
        cat("", file=filename.gz)
    }

    for (imod in 1:length(hmm$info$ID)) {
        ID <- hmm$info$ID[imod]
        cond <- hmm$info$condition[imod]
        
        ## Select only segments with peaks
        peaks <- insertchr(hmm$peaks[[ID]])
        peaks.df <- as.data.frame(peaks)
        peaks.df$peakNumber <- paste0('peak_', 1:nrow(peaks.df))
        
        # Data.frame for write.table
        df <- peaks.df[,c('chromosome','start','end','peakNumber','peakScores','strand')]
        df$strand <- sub('\\*', '.', df$strand)
        
        # Make score integer
        df$peakScores <- round(df$peakScores*1000)
        
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$start
        df$thickEnd <- df$end
        # Colors
        RGB <- t(grDevices::col2rgb(getStateColors('modified')))
        df$itemRgb <- apply(RGB,1,paste,collapse=",")
        
        ## Write to file
        if (separate.files) {
            filename.sep <- paste0(sub('.bed.gz$', '', filename), '_', ID, '.bed.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('Writing to file ',filename, ' ...')
        }
        if (header) {
            if (is.null(trackname)) {
                trackname.string <- paste0("peaks for ", ID)
            } else {
                trackname.string <- paste0("peaks for ", ID, ", ", trackname)
            }
            priority <- 52 + 4*imod
            cat(paste0('track name="', trackname.string, '" description="', trackname.string, '" visibility=1 itemRgb=On priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        if (separate.files) {
            close(filename.gz)
        }
        stopTimedMessage(ptm)
    }
    if (!separate.files) {
        close(filename.gz)
    }
}

#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportCombinedMultivariateCombinations <- function(hmm, filename, exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=TRUE) {

    ## Intercept user input
    if (is.data.frame(exclude.states)) {
        exclude.states <- exclude.states$combination
    }
    if (is.data.frame(include.states)) {
        include.states <- include.states$combination
    }

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.combined.multivariate.hmm)[[1]]

    ## Write first line to file
    if (!separate.files) {
        filename <- paste0(filename,".bed.gz")
        filename.gz <- gzfile(filename, 'w')
        ptm <- startTimedMessage('Writing to file ',filename, ' ...')
        cat("", file=filename.gz)
    }

    ## Export the combinatorial states from individual segments
    conditions <- names(hmm$segments.per.condition)
    for (cond in conditions) {
        segments <- hmm$segments.per.condition[[cond]]

        # Exclude and include states
        combstates <- levels(segments$combination)
        if (is.null(include.states)) {
            combstates2use <- setdiff(combstates, exclude.states)
        } else {
            combstates2use <- intersect(combstates, include.states)
        }
        segments <- segments[mcols(segments)$combination %in% combstates2use]
        segments <- insertchr(segments)
        segments.df <- as.data.frame(segments)

        # Make header
        if (is.null(trackname)) {
            trackname.cond <- paste0('condition: ', cond, ', combinations')
        } else {
            trackname.cond <- paste0('condition: ', cond, ', combinations, ', trackname)
        }
        # Make filenames
        if (separate.files) {
            filename.cond <- paste0(filename, '_', cond, '.bed.gz')
            filename.gz <- gzfile(filename.cond, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.cond, ' ...')
            cat("", file=filename.gz)
        }

        ## Export combinatorial states
        # Generate the colors for each combinatorial state
        colors <- getDistinctColors(length(levels(segments.df$combination)))
        RGBs <- t(grDevices::col2rgb(colors))
        RGBlist <- list()
        for (i1 in 1:3) {
            RGBlist[[i1]] <- RGBs[,i1]
        }
        RGBlist$sep <- ','
        RGBs <- do.call(paste, RGBlist)
        itemRgb <- RGBs[as.integer(factor(segments.df$combination, levels=sort(levels(segments.df$combination))))]

        # Write to file
        numsegments <- nrow(segments.df)
        if (is.null(segments.df$differential.score)) {
            segments.df$differential.score <- 0
        }
        names(segments.df)[names(segments.df)=='differential.score'] <- 'score'
        df <- cbind(segments.df[,c('chromosome','start','end','combination','score')], strand=rep(".",numsegments), thickStart=segments.df$start, thickEnd=segments.df$end, itemRgb=itemRgb)
        # Make score integer
        df$score <- round(df$score*1000)
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$thickStart - 1
        if (header) {
            cat(paste0('track name="',trackname.cond, '" description="', trackname.cond, '" visibility=1 itemRgb=On priority=49\n'), file=filename.gz, append=TRUE)
        }
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

        if (separate.files) {
            close(filename.gz)
            stopTimedMessage(ptm)
        }
    }
    if (!separate.files) {
        close(filename.gz)
        stopTimedMessage(ptm)
    }

}


#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportCombinedMultivariateCounts <- function(hmm, filename, header=TRUE, separate.files=TRUE, trackname=NULL) {

    ## Load models
    hmm <- loadHmmsFromFiles(hmm, check.class=class.combined.multivariate.hmm)[[1]]

    ## RPKM normalization
    hmm$bins$counts <- rpkm.matrix(hmm$bins$counts, binsize=mean(width(hmm$bins)))
    
    
    ## Variables
    filename <- paste0(filename,".wig.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }
    nummod <- length(hmm$info$ID)
    readcol <- paste(grDevices::col2rgb(getStateColors('counts')), collapse=',')

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        ID <- hmm$info$ID[imod]
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', ID, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('  Writing track ',imod,' / ',nummod, ' ...')
        }
        priority <- 50 + 4*imod
        binsize <- width(hmm$bins[1])
        if (header) {
            cat(paste0('track type=wiggle_0 name="read count for ',ID,'" description="read count for ',ID,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        for (chrom in seqlevels(hmm$bins)) {
            if (!grepl('chr', chrom)) {
                chromr <- sub(pattern='^', replacement='chr', chrom)
            } else {
                chromr <- chrom
            }
            cat(paste0("fixedStep chrom=",chromr," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(hmm$bins[seqnames(hmm$bins)==chrom]$counts[,imod], file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
        }
        stopTimedMessage(ptm)
    }
    if (!separate.files) {
        close(filename.gz)
    }
}


#====================================================
# Export regions from GRanges
#====================================================
#' Export genome browser viewable files
#'
#' Export GRanges as genome browser viewable file
#'
#' Export regions from \code{\link{GRanges}} as a file which can be uploaded into a genome browser. Regions are exported in BED format (.bed.gz).
#'
#' @author Aaron Taudt
#' @param gr A \code{\link{GRanges}} object.
#' @param trackname The name that will be used as track name and description in the header.
#' @param filename The name of the file that will be written. The ending ".bed.gz". Any existing file will be overwritten.
#' @param namecol A character specifying the column that is used as name-column.
#' @param scorecol A character specifying the column that is used as score-column. The score should contain integers in the interval [0,1000] for compatibility with the UCSC genome browser convention.
#' @param colorcol A character specifying the column that is used for coloring the track. There will be one color for each unique element in \code{colorcol}.
#' @param colors A character vector with the colors that are used for the unique elements in \code{colorcol}.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param append Whether or not to append to an existing file.
#' @return \code{NULL}
#' @seealso \code{\link{exportPeaks}}, \code{\link{exportCounts}}, \code{\link{exportCombinations}}
#' @importFrom utils write.table
#' @export
#' @examples 
#'### Export regions with read counts above 20 ###
#'# Get an example BAM file with ChIP-seq reads
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'# Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(file, assembly=rn4_chrominfo, binsizes=1000,
#'                   chromosomes='chr12')
#'plotHistogram(binned)
#'# Export regions with read count above 20
#'exportGRangesAsBedFile(binned[binned$counts > 20], filename=tempfile(),
#'              trackname='read counts above 20')
#'
exportGRangesAsBedFile <- function(gr, trackname, filename, namecol='combination', scorecol='score', colorcol=NULL, colors=NULL, header=TRUE, append=FALSE) {

    if (length(gr)==0) {
        warning("Supplied GRanges object contains no ranges.")
        return()
    }

    ## Transform to GRanges
    gr <- insertchr(gr)

    # Variables
    filename <- paste0(filename,".bed.gz")
    if (append) {
        filename.gz <- gzfile(filename, 'a')
    } else {
        filename.gz <- gzfile(filename, 'w')
    }

    # Write first line to file
    if (append) {
        ptm <- startTimedMessage('Appending to file ',filename, ' ...')
    } else {
        ptm <- startTimedMessage('Writing to file ',filename, ' ...')
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    if (header) {
        cat(paste0("track name=\"",trackname,"\" description=\"",trackname,"\" visibility=1 itemRgb=On\n"), file=filename.gz, append=TRUE)
    }
    if (! scorecol %in% names(mcols(gr))) {
        gr$score <- 0
    } else {
        gr$score <- mcols(gr)[,scorecol]
				# Check if scorecolumn follows the UCSC convention
				if (min(gr$score) < 0 | max(gr$score) > 1000 | all(gr$score<=1.2) | !is.integer(gr$score)) {
						warning("Column '", scorecol, "' should contain integer values between 0 and 1000 for compatibility with the UCSC convention.")
				}
    }
    regions <- as.data.frame(gr)[c('chromosome','start','end','score','strand')]
    regions$strand <- sub("\\*", ".", regions$strand)
    if (! namecol %in% names(mcols(gr))) {
        regions$name <- paste0('region_', 1:nrow(regions))
    } else {
        regions$name <- mcols(gr)[,namecol]
    }
    regions <- regions[c('chromosome','start','end','name','score','strand')]
    # Convert from 1-based closed to 0-based half open
    regions$start <- regions$start - 1
    df <- regions
    
    if (!is.null(colorcol)) {
        df <- cbind(regions, thickStart=regions$start, thickEnd=regions$end)
        # Generate the colors for each element in 'namecol'
        if (colorcol %in% names(mcols(gr))) {
            if (is.null(colors)) {
                colors <- getDistinctColors(length(unique(df$name)))
            }
            RGBs <- t(grDevices::col2rgb(colors))
            RGBlist <- list()
            for (i1 in 1:3) {
                RGBlist[[i1]] <- RGBs[,i1]
            }
            RGBlist$sep <- ','
            RGBs <- do.call(paste, RGBlist)
            itemRgb <- RGBs[as.integer(factor(df$name))]
            df$itemRgb <- itemRgb
        }
    }
    
    if (nrow(df) == 0) {
        warning('No regions in input')
    } else {
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    }

    close(filename.gz)
    stopTimedMessage(ptm)

}
