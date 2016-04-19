#=====================================================
# Export univariate HMMs
#=====================================================
#' Export genome browser viewable files
#'
#' Export univariate peak-calls and read counts as genome browser viewable file
#'
#' Export \code{\link{uniHMM}} objects as files which can be uploaded into a genome browser. Peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format (.wig.gz).
#'
#' @author Aaron Taudt
#' @param hmm.list A \code{list()} of \code{\link{uniHMM}} objects or vector of files that contain such objects.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for peak-calls or ".wig.gz" for read counts. Any existing file will be overwritten.
#' @param what A character vector specifying what will be exported. Supported are \code{c('peaks', 'counts')}.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param separate.files A logical indicating whether or not to produce separate files for each hmm in \code{hmm.list}.
#' @param orderByScore Logical indicating whether or not to sort entries in BED file by score.
#' @return \code{NULL}
#' @seealso \code{\link{exportBinnedData}}, \code{\link{exportMultivariate}}
#' @export
#' @examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(bedfile, assembly=rn4_chrominfo, binsize=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned, ID='example_H3K27me3', max.time=60, eps=1)
#'## Export
#'exportUnivariates(list(hmm), filename=tempfile(), what=c('peaks','counts'))
#'
exportUnivariates <- function(hmm.list, filename, what=c('peaks', 'counts'), header=TRUE, separate.files=FALSE, orderByScore=TRUE) {
    if ('peaks' %in% what) {
        exportUnivariatePeaks(hmm.list, filename=paste0(filename, '_peaks'), header=header, separate.files=separate.files, orderByScore=orderByScore)
    }
    if ('counts' %in% what) {
        exportUnivariateReadCounts(hmm.list, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files)
    }
}

#----------------------------------------------------
# Export peak-calls from univariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportUnivariatePeaks <- function(hmm.list, filename="chromstaR_univariatePeakCalls", header=TRUE, separate.files=FALSE, orderByScore=TRUE) {

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

    ## Load models
    hmm.list <- loadHmmsFromFiles(hmm.list, check.class=class.univariate.hmm)

    ## Transform to GRanges
    hmm.grl <- lapply(hmm.list, '[[', 'segments')
    hmm.grl <- lapply(hmm.grl, insertchr)

    # Variables
    nummod <- length(hmm.list)
    filename <- paste0(filename,".bed.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }

    # Generate the colors
    colors <- getStateColors(levels(hmm.grl[[1]]$state))
    RGBs <- t(grDevices::col2rgb(colors))
    RGBs <- apply(RGBs,1,paste,collapse=",")

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        message('  Writing track ',imod,' / ',nummod)
        hmm <- hmm.list[[imod]]
        if (separate.files) {
            filename.sep <- paste0(sub('.bed.gz$', '', filename), '_', hmm$ID, '.bed.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            message('Writing to file ',filename.sep)
            cat("", file=filename.gz)
        }
        hmm.gr <- hmm.grl[[imod]]
        priority <- 51 + 4*imod
        if (header) {
            cat(paste0("track name=\"univariate calls for ",hmm$ID,"\" description=\"univariate calls for ",hmm$ID,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
        }
        if (is.null(hmm.gr$score)) {
            hmm.gr$score <- 0
        }
        collapsed.calls <- as.data.frame(hmm.gr)[c('chromosome','start','end','state','score')]
        # Make score integer
        collapsed.calls$score <- round(collapsed.calls$score*1000)
        itemRgb <- RGBs[as.character(collapsed.calls$state)]
        numsegments <- nrow(collapsed.calls)
        df <- cbind(collapsed.calls, strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
        # Reorder
        if (orderByScore) {
            df <- df[rev(order(df$score)),]
        }
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$thickStart - 1
        df <- df[df$state=='modified',]
        if (nrow(df) == 0) {
            warning('hmm ',imod,' does not contain any \'modified\' calls')
        } else {
            utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
        }
    }
    if (!separate.files) {
        close(filename.gz)
    }

}


#----------------------------------------------------
# Export read counts from univariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportUnivariateReadCounts <- function(hmm.list, filename="chromstaR_univariateReadCounts", header=TRUE, separate.files=FALSE) {

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

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
        message('  Writing track ',imod,' / ',nummod)
        hmm <- hmm.list[[imod]]
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', hmm$ID, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            message('Writing to file ',filename.sep)
            cat("", file=filename.gz)
        }
        hmm.gr <- hmm.grl[[imod]]
        priority <- 50 + 4*imod
        binsize <- width(hmm.gr[1])
        if (header) {
            cat(paste0('track type=wiggle_0 name="read count for ',hmm$ID,'" description="read count for ',hmm$ID,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        for (chrom in unique(hmm.gr$chromosome)) {
            cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
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
# Export multivariate HMMs
#=====================================================
#' Export genome browser viewable files
#'
#' Export multivariate calls and read counts as genome browser viewable file
#'
#' Export \code{\link{uniHMM}} objects as files which can be uploaded into a genome browser. Combinatorial states and peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format (.wig.gz).
#'
#' @author Aaron Taudt
#' @param multi.hmm A \code{\link{multiHMM}} object or file that contains such an object.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for combinatorial states and peak-calls or ".wig.gz" for read counts. Any existing file will be overwritten.
#' @param what A character vector specifying what will be exported. Supported are \code{c('combinations', 'peaks', 'counts')}.
#' @param exclude.states A character vector with combinatorial states that will be excluded from export.
#' @param include.states A character vector with combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
#' @param trackname Name that will be used in the "track name" field of the BED file.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param separate.files A logical indicating whether or not to produce separate files for peaks if \code{what} contains 'peaks' or 'counts'.
#' @param orderByScore Logical indicating whether or not to sort entries in BED file by score.
#' @return \code{NULL}
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportBinnedData}}
#' @export
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-mark_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'## Export peak calls and combinatorial states
#'exportMultivariate(model, filename=tempfile(), what=c('peaks','combinations'))
#'
exportMultivariate <- function(multi.hmm, filename, what=c('combinations', 'peaks', 'counts'), exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=FALSE, orderByScore=TRUE) {
    if ('combinations' %in% what) {
        exportMultivariateCalls(multi.hmm, filename=paste0(filename, '_combinations'), separate.tracks=FALSE, exclude.states, include.states, trackname=trackname, header=header, orderByScore=orderByScore)
    }
    if ('peaks' %in% what) {
        exportMultivariateCalls(multi.hmm, filename=paste0(filename, '_peaks'), separate.tracks=TRUE, exclude.states, include.states, header=header, separate.files=separate.files, orderByScore=orderByScore)
    }
    if ('counts' %in% what) {
        exportMultivariateReadCounts(multi.hmm, filename=paste0(filename, '_counts'), header=header, separate.files=separate.files)
    }
}

#----------------------------------------------------
# Export combinatorial states or peak-calls from multivariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportMultivariateCalls <- function(multi.hmm, filename="chromstaR_multivariateCalls", separate.tracks=TRUE, exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=FALSE, orderByScore=TRUE) {

    ## Intercept user input
    if (class(multi.hmm)!=class.multivariate.hmm) {
        multi.hmm <- get(load(multi.hmm))
        if (class(multi.hmm)!=class.multivariate.hmm) {
            stop("argument 'multi.hmm' expects a multivariate hmm or a file which contains a multivariate hmm")
        }
    }

    if (is.data.frame(exclude.states)) {
        exclude.states <- exclude.states$combination
    }
    if (is.data.frame(include.states)) {
        include.states <- include.states$combination
    }

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

    ## Variables
    combstates <- levels(multi.hmm$bins$combination)
    numstates <- length(combstates)
    if (is.null(include.states)) {
        combstates2use <- setdiff(combstates, exclude.states)
    } else {
        combstates2use <- intersect(combstates, include.states)
    }
    numstates <- length(combstates2use)
    nummod <- length(multi.hmm$IDs)

    ## Insert chromosome if missing
    segments.df <- as.data.frame(insertchr(multi.hmm$segments))

    # Select only desired states
    segments.df <- segments.df[segments.df$combination %in% combstates2use,]
    if (nrow(segments.df) == 0) {
        warning("No regions to export!")
        return()
    }

    # Write first line to file
    filename <- paste0(filename,".bed.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }

    if (separate.tracks) {
        ## Export peak calls
        bin <- dec2bin(segments.df$state, ndigits=length(multi.hmm$IDs))
        colnames(bin) <- multi.hmm$IDs
        for (imod in 1:nummod) {
            ID <- multi.hmm$IDs[imod]
            if (separate.files) {
                filename.sep <- paste0(sub('.bed.gz$', '', filename), '_', ID, '.bed.gz')
                filename.gz <- gzfile(filename.sep, 'w')
                ptm <- startTimedMessage('  Writing to file ',filename.sep, ' ...')
                cat("", file=filename.gz)
            } else {
                ptm <- startTimedMessage('  Writing track ',imod,' / ',nummod, ' ...')
            }
            numsegments <- length(which(bin[,imod]))
            priority <- 52 + 4*imod
            mask <- bin[,imod]
            if (length(grep('^mean.posteriors', names(segments.df)))==ncol(bin)) {
                segments.df$score <- segments.df[,grep('^mean.posteriors', names(segments.df))[imod]]
            } else {
                segments.df$score <- 0
            }
            if (is.null(segments.df$combination)) {
                df <- cbind(segments.df[mask,c('chromosome','start','end','state','score')], strand=rep(".",numsegments))
            } else {
                df <- cbind(segments.df[mask,c('chromosome','start','end','combination','score')], strand=rep(".",numsegments))
            }
            # Make score integer
            df$score <- round(df$score*1000)
            # Reorder
            if (orderByScore) {
                df <- df[rev(order(df$score)),]
            }
            # Convert from 1-based closed to 0-based half open
            df$start <- df$start - 1
            df$thickStart <- df$start
            df$thickEnd <- df$end
            RGB <- t(grDevices::col2rgb(getStateColors('modified')))
            RGB <- apply(RGB,1,paste,collapse=",")
            df$itemRgb <- rep(RGB, numsegments)
            if (header) {
                cat(paste0("track name=\"multivariate peaks for ",colnames(bin)[imod],"\" description=\"multivariate peaks for ",colnames(bin)[imod],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
            }
            utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
            if (separate.files) {
                close(filename.gz)
            }
            stopTimedMessage(ptm)
        }
    } else {
        ## Export combinatorial states
        # Generate the colors for each combinatorial state
        colors <- getDistinctColors(length(levels(segments.df$combination)))
        RGBs <- t(grDevices::col2rgb(colors))
        RGBs <- apply(RGBs,1,paste,collapse=",")
        itemRgb <- RGBs[as.integer(factor(segments.df$combination, levels=sort(levels(segments.df$combination))))]

        # Write to file
        numsegments <- nrow(segments.df)
        if (is.null(segments.df$score)) {
            segments.df$score <- 0
        }
        if (is.null(segments.df$combination)) {
            df <- cbind(segments.df[,c('chromosome','start','end','state','score')], strand=rep(".",numsegments), thickStart=segments.df$start, thickEnd=segments.df$end, itemRgb=itemRgb)
        } else {
            df <- cbind(segments.df[,c('chromosome','start','end','combination','score')], strand=rep(".",numsegments), thickStart=segments.df$start, thickEnd=segments.df$end, itemRgb=itemRgb)
        }
        # Make score integer
        df$score <- round(df$score*1000)
        # Reorder
        if (orderByScore) {
            df <- df[rev(order(df$score)),]
        }
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$thickStart - 1
        if (header) {
            if (is.null(trackname)) {
                cat(paste0('track name="combinations" description="combinations" visibility=1 itemRgb=On priority=49\n'), file=filename.gz, append=TRUE)
            } else {
                cat(paste0('track name="',trackname,'" description="',trackname,'" visibility=1 itemRgb=On priority=49\n'), file=filename.gz, append=TRUE)
            }
        }
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    }
    if (!separate.files) {
        close(filename.gz)
    }

}


#----------------------------------------------------
# Export read counts from multivariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportMultivariateReadCounts <- function(multi.hmm, filename="chromstaR_multivariateReadCounts", header=TRUE, separate.files=FALSE) {

    if (class(multi.hmm)!=class.multivariate.hmm) {
        multi.hmm <- get(load(multi.hmm))
        if (class(multi.hmm)!=class.multivariate.hmm) {
            stop("argument 'multi.hmm' expects a multivariate hmm or a file which contains a multivariate hmm")
        }
    }

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

    ## Variables
    filename <- paste0(filename,".wig.gz")
    if (!separate.files) {
        filename.gz <- gzfile(filename, 'w')
    }
    nummod <- length(multi.hmm$IDs)
    readcol <- paste(grDevices::col2rgb(getStateColors('counts')), collapse=',')

    # Write first line to file
    if (!separate.files) {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    for (imod in 1:nummod) {
        ID <- multi.hmm$IDs[imod]
        if (separate.files) {
            filename.sep <- paste0(sub('.wig.gz$', '', filename), '_', ID, '.wig.gz')
            filename.gz <- gzfile(filename.sep, 'w')
            ptm <- startTimedMessage('  Writing to file ',filename.sep, ' ...')
            cat("", file=filename.gz)
        } else {
            ptm <- startTimedMessage('  Writing track ',imod,' / ',nummod, ' ...')
        }
        priority <- 50 + 4*imod
        binsize <- width(multi.hmm$bins[1])
        if (header) {
            cat(paste0('track type=wiggle_0 name="read count for ',ID,'" description="read count for ',ID,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        for (chrom in seqlevels(multi.hmm$bins)) {
            if (!grepl('chr', chrom)) {
                chromr <- sub(pattern='^', replacement='chr', chrom)
            } else {
                chromr <- chrom
            }
            cat(paste0("fixedStep chrom=",chromr," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(multi.hmm$bins[seqnames(multi.hmm$bins)==chrom]$counts[,imod], file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
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
# Export binned data
#=====================================================
#' Export genome browser viewable files
#'
#' Export read counts as genome browser viewable file
#'
#' Export read counts from \code{\link{binned.data}} as a file which can be uploaded into a genome browser. Read counts are exported in WIGGLE format (.wig.gz).
#'
#' @author Aaron Taudt
#' @param binned.data.list A \code{list()} of \code{\link{binned.data}} objects or vector of files that contain such objects.
#' @param filename The name of the file that will be written. The ending ".wig.gz" for read counts will be appended. Any existing file will be overwritten.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param separate.files A logical indicating whether or not to produce separate files for each object in \code{binned.data.list}.
#' @return \code{NULL}
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportMultivariate}}
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
#' @examples
#'## Get example BED-files with ChIP-seq reads for H3K36me3
#' # in 7 different brain tissues (chr22)
#'bedfiles <- list.files(system.file(file.path("extdata","brain"),
#'                       package="chromstaR"), full=TRUE)
#'## Bin the data into bin size 1000bp
#'binned.data.list <- list()
#'for (bedfile in bedfiles) {
#'  binned.data.list[[bedfile]] <- bed2binned(bedfile, assembly='hg19', binsize=1000,
#'                                            save.as.RData=FALSE)
#'}
#'## Export the binned read counts
#'\donttest{exportBinnedData(binned.data.list, filename='chromstaR-example_binned.data')}
#' @export
exportBinnedData <- function(binned.data.list, filename="chromstaR_ReadCounts", header=TRUE, separate.files=FALSE) {

    ## Load data
    if (is.character(binned.data.list)) {
        ptm <- startTimedMessage('Loading binned.data from files ...')
        binfiles <- binned.data.list
        binned.data.list <- list()
        for (binfile in binfiles) {
            binned.data.list[[binfile]] <- get(load(binfile))
        }
        stopTimedMessage(ptm)
    }

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

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
            cat(paste0('track type=wiggle_0 name="read count for ',name,'" description="read count for ',name,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
        }
        # Write read data
        for (chrom in unique(b$chromosome)) {
            cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
            utils::write.table(mcols(b[b$chromosome==chrom])$counts, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
        if (separate.files) {
            close(filename.gz)
        }
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
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param orderByScore If set to \code{TRUE}, the lines in the output will be ordered by column 'score' in descending order.
#' @param append Whether or not to append to an existing file.
#' @return \code{NULL}
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportMultivariate}}
#' @importFrom utils write.table
#' @export
#' @examples 
#'### Export regions with read counts above 20 ###
#'# Get an example BED file with ChIP-seq reads
#'bedfile <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'# Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(bedfile, assembly=rn4_chrominfo, binsize=1000,
#'                   chromosomes='chr12')
#'plot(binned)
#'# Export regions with read count above 20
#'exportGRanges(binned[binned$counts > 20], filename=tempfile(),
#'              trackname='read counts above 20')
#'
exportGRanges <- function(gr, trackname, filename="chromstaR_GRanges_regions", header=TRUE, orderByScore=FALSE, append=FALSE) {

    if (length(gr)==0) {
        warning("Supplied GRanges object contains no ranges.")
        return()
    }

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
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
        message('appending to file ',filename)
    } else {
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }
    
    ### Write every model to file ###
    if (header) {
        cat(paste0("track name=\"",trackname,"\" description=\"",trackname,"\" visibility=1 itemRgb=Off\n"), file=filename.gz, append=TRUE)
    }
    if (is.null(gr$score)) {
        gr$score <- 0
    } else {
        if (orderByScore) {
            gr <- gr[rev(order(gr$score))]
        }
    }
    regions <- as.data.frame(gr)[c('chromosome','start','end','score')]
    regions$name <- paste0('region_', 1:nrow(regions))
    regions <- regions[c('chromosome','start','end','name','score')]
    # Make score integer
    regions$score <- round(regions$score*1000)
    numsegments <- nrow(regions)
    df <- cbind(regions, strand=rep(".",numsegments), thickStart=regions$start, thickEnd=regions$end)
    # Convert from 1-based closed to 0-based half open
    df$start <- df$start - 1
    df$thickStart <- df$thickStart - 1
    if (nrow(df) == 0) {
        warning('No regions in input')
    } else {
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    }

    close(filename.gz)

}


#=====================================================
# Export combined multivariate HMMs
#=====================================================
#' Export genome browser viewable files
#'
#' Export multivariate calls as genome browser viewable file
#'
#' Export \code{\link{combinedMultiHMM}} objects as files which can be uploaded into a genome browser. Combinatorial states are exported in BED format (.bed.gz).
#'
#' @author Aaron Taudt
#' @param hmm A \code{\link{combinedMultiHMM}} object or file that contains such an object.
#' @param filename The name of the file that will be written. The ending ".bed.gz" for combinatorial states will be appended. Any existing file will be overwritten.
#' @param what A character vector specifying what will be exported. Supported are \code{c('combinations')}.
#' @param exclude.states A vector of combinatorial states that will be excluded from export.
#' @param include.states A vector of combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
#' @param trackname Name that will be used in the "track name" field of the BED file.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @param separate.files A logical indicating whether or not to produce separate files for each condition.
#' @return \code{NULL}
#' @export
exportCombinedMultivariate <- function(hmm, filename, what=c('combinations'), exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=FALSE) {
    if ('combinations' %in% what) {
        exportCombinedMultivariateCalls(hmm, filename=filename, exclude.states=exclude.states, include.states=include.states, trackname=trackname, header=header, separate.files=separate.files)
    }
}

#----------------------------------------------------
# Export combinatorial states or peak-calls from multivariate HMMs
#----------------------------------------------------
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
exportCombinedMultivariateCalls <- function(hmm, filename="chromstaR_combinedMultivariateCalls", exclude.states='[]', include.states=NULL, trackname=NULL, header=TRUE, separate.files=FALSE) {

    if (class(hmm)!=class.combined.multivariate.hmm) {
        hmm <- get(load(hmm))
        if (class(hmm)!=class.combined.multivariate.hmm) {
            stop("argument 'hmm' expects a combined multivariate hmm or a file which contains an object")
        }
    }

    ## Function definitions
    insertchr <- function(hmm.gr) {
        # Change chromosome names from '1' to 'chr1' if necessary
        mask <- which(!grepl('chr', seqnames(hmm.gr)))
        mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
        mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
        mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
        return(hmm.gr)
    }

    ## Write first line to file
    if (!separate.files) {
        filename <- paste0(filename,".bed.gz")
        filename.gz <- gzfile(filename, 'w')
        message('Writing to file ',filename)
        cat("", file=filename.gz)
    }

    ## Export the combinatorial states from individual segments
    conditions <- names(hmm$segments.separate)
    for (cond in conditions) {
        segments <- hmm$segments.separate[[cond]]

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
            trackname.cond <- paste0('condition: ', cond, ', combinatorial states')
        } else {
            trackname.cond <- paste0('condition: ', cond, ', ', trackname)
        }
        # Make filenames
        if (separate.files) {
            filename.cond <- paste0(filename, '_', cond, '.bed.gz')
            filename.gz <- gzfile(filename.cond, 'w')
            message('Writing to file ',filename.cond)
            cat("", file=filename.gz)
        }

        ## Export combinatorial states
        # Generate the colors for each combinatorial state
        colors <- getDistinctColors(length(levels(segments.df$combination)))
        RGBs <- t(grDevices::col2rgb(colors))
        RGBs <- apply(RGBs,1,paste,collapse=",")
        itemRgb <- RGBs[as.integer(factor(segments.df$combination, levels=sort(levels(segments.df$combination))))]

        # Write to file
        numsegments <- nrow(segments.df)
        if (is.null(segments.df$score)) {
            segments.df$score <- 0
        }
        df <- cbind(segments.df[,c('chromosome','start','end','combination','score')], strand=rep(".",numsegments), thickStart=segments.df$start, thickEnd=segments.df$end, itemRgb=itemRgb)
        # Make score integer
        df$score <- round(df$score*1000)
        # Reorder
#         if (orderByScore) {
#             df <- df[rev(order(df$score)),]
#         }
        # Convert from 1-based closed to 0-based half open
        df$start <- df$start - 1
        df$thickStart <- df$thickStart - 1
        if (header) {
            cat(paste0('track name="',trackname.cond, '" description="', trackname.cond, '" visibility=1 itemRgb=On priority=49\n'), file=filename.gz, append=TRUE)
        }
        utils::write.table(format(df, scientific=FALSE, trim=TRUE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

        if (separate.files) {
            close(filename.gz)
        }
    }
    if (!separate.files) {
        close(filename.gz)
    }

}
