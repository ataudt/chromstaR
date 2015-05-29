#=====================================================
# Export univariate HMMs
#=====================================================
#' Export genome browser viewable files
#'
#' Export univariate peak-calls and read counts as genome browser viewable file
#'
#' Export \code{\link{chromstaR_univariateHMM}} objects as files which can be uploaded into a genome browser. Peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format (.wiggle.gz).
#'
#' @author Aaron Taudt
#' @param hmm.list A list of \code{\link{chromstaR_univariateHMM}} objects or files that contain such objects.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for peak-calls or ".wiggle.gz" for read counts. Any existing file will be overwritten.
#' @param what A character vector specifying what will be exported. Supported are \code{c('peaks', 'reads')}.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @seealso \code{\link{exportBinnedData}}, \code{\link{exportMultivariate}}
#' @examples
#'### Univariate peak calling ###
#'## Get example BED-files with ChIP-seq reads for H3K36me3
#' # in 7 different brain tissues (chr22)
#'bedfiles <- list.files(system.file(file.path("extdata","brain"),
#'                       package="chromstaR"), full=TRUE)
#'## Bin the data into bin size 1000bp and build the univariate Hidden Markov Model (HMM)
#'binned.data.list <- list()
#'uni.HMMs <- list()
#'for (bedfile in bedfiles) {
#'  binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000, save.as.RData=FALSE)
#'  binned.data.list[[bedfile]] <- binned.data
#'  uni.HMMs[[bedfile]] <- callPeaksUnivariate(binned.data, ID=basename(bedfile),
#'                                             max.time=30, eps=0.01)
#'}
#'## Export the binned read counts and peaks
#'\donttest{exportUnivariates(uni.HMMs, filename='chromstaR-example_univariate.HMMs', what=c('reads','peaks'))}
#' @export
exportUnivariates <- function(hmm.list, filename, what=c('peaks', 'reads'), header=TRUE) {
	if ('peaks' %in% what) {
		exportUnivariatePeaks(hmm.list, filename, header=header)
	}
	if ('reads' %in% what) {
		exportUnivariateReadCounts(hmm.list, filename, header=header)
	}
}

#----------------------------------------------------
# Export peak-calls from univariate HMMs
#----------------------------------------------------
exportUnivariatePeaks <- function(hmm.list, filename="chromstaR_univariatePeakCalls", header=TRUE) {

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
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges
	hmm.grl <- lapply(hmm.list, '[[', 'segments')
	hmm.grl <- lapply(hmm.grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	filename <- paste0(filename,".bed.gz")
	filename.gz <- gzfile(filename, 'w')

	# Generate the colors
	colors <- state.colors[levels(hmm.grl[[1]]$state)]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")

	# Write first line to file
	message('writing to file ',filename)
# 	cat("browser hide all\n", file=filename.gz)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 51 + 4*imod
		if (header) {
			cat(paste0("track name=\"univariate calls for ",hmm$ID,"\" description=\"univariate calls for ",hmm$ID,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
		}
		collapsed.calls <- as.data.frame(hmm.gr)[c('chromosome','start','end','state','mean.posterior.modified')]
		itemRgb <- RGBs[as.character(collapsed.calls$state)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		df$thickStart <- df$thickStart - 1
		df <- df[df$state=='modified',]
		if (nrow(df) == 0) {
			warning('hmm ',imod,' does not contain any \'modified\' calls')
		} else {
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}
	close(filename.gz)
	message('')

}


#----------------------------------------------------
# Export read counts from univariate HMMs
#----------------------------------------------------
exportUnivariateReadCounts <- function(hmm.list, filename="chromstaR_univariateReadCounts", header=TRUE) {

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
	hmm.list <- loadHmmsFromFiles(hmm.list)

	## Transform to GRanges
	grl <- lapply(hmm.list, '[[', 'bins')
	hmm.grl <- lapply(grl, insertchr)

	# Variables
	nummod <- length(hmm.list)
	filename <- paste0(filename,".wiggle.gz")
	filename.gz <- gzfile(filename, 'w')
	readcol <- paste(col2rgb(state.colors['reads']), collapse=',')

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 50 + 4*imod
		binsize <- width(hmm.gr[1])
		if (header) {
			cat(paste0('track type=wiggle_0 name="read count for ',hmm$ID,'" description="read count for ',hmm$ID,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
		}
		# Write read data
		for (chrom in unique(hmm.gr$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$reads, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
}


#=====================================================
# Export multivariate HMMs
#=====================================================
#' Export genome browser viewable files
#'
#' Export multivariate calls and read counts as genome browser viewable file
#'
#' Export \code{\link{chromstaR_univariateHMM}} objects as files which can be uploaded into a genome browser. Combinatorial states and peak-calls are exported in BED format (.bed.gz) and read counts are exported in WIGGLE format (.wiggle.gz).
#'
#' @author Aaron Taudt
#' @param multi.hmm A \code{\link{chromstaR_multivariateHMM}} object or file that contains such an object.
#' @param filename The name of the file that will be written. The appropriate ending will be appended, either ".bed.gz" for combinatorial states and peak-calls or ".wiggle.gz" for read counts. Any existing file will be overwritten.
#' @param what A character vector specifying what will be exported. Supported are \code{c('combstates', 'peaks', 'reads')}.
#' @param exclude.states A vector of combinatorial states that will be excluded from export.
#' @param include.states A vector of combinatorial states that will be exported. If specified, \code{exclude.states} is ignored.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportBinnedData}}
#' @examples
#'### Multivariate peak calling ###
#'## Load example multivariate Hidden Markov Model
#'data(example.multi.HMM)
#'## Export the binned read counts and peaks for each track
#'\donttest{
#'exportMultivariate(example.multi.HMM, exclude.states=c(0),
#'                   filename='chromstaR_example.multi.HMM',
#'                   what=c('reads','peaks'))
#'}
#' @export
exportMultivariate <- function(multi.hmm, filename, what=c('combstates', 'peaks', 'reads'), exclude.states=0, include.states=NULL, header=TRUE) {
	if ('combstates' %in% what) {
		exportMultivariateCalls(multi.hmm, filename, separate.tracks=FALSE, exclude.states, include.states, header=header)
	}
	if ('peaks' %in% what) {
		exportMultivariateCalls(multi.hmm, filename, separate.tracks=TRUE, exclude.states, include.states, header=header)
	}
	if ('reads' %in% what) {
		exportMultivariateReadCounts(multi.hmm, filename, header=header)
	}
}

#----------------------------------------------------
# Export combinatorial states or peak-calls from multivariate HMMs
#----------------------------------------------------
exportMultivariateCalls <- function(multi.hmm, filename="chromstaR_multivariateCalls", separate.tracks=TRUE, exclude.states=0, include.states=NULL, header=TRUE) {

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
	filename <- paste0(filename,".bed.gz")
	filename.gz <- gzfile(filename, 'w')
	combstates <- levels(multi.hmm$bins$state)
	numstates <- length(combstates)
	if (is.null(include.states)) {
		combstates2use <- setdiff(combstates, exclude.states)
	} else {
		combstates2use <- intersect(combstates, include.states)
	}
	numstates <- length(combstates2use)

	## Collapse the calls
	if (is.null(multi.hmm$segments$combination)) {
		collapsed.calls <- as.data.frame(insertchr(multi.hmm$segments))[,c('chromosome','start','end','state')]
	} else {
		collapsed.calls <- as.data.frame(insertchr(multi.hmm$segments))[,c('chromosome','start','end','state','combination')]
	}
	# Select only desired states
	mask <- rep(FALSE,nrow(collapsed.calls))
	for (istate in combstates2use) {
		mask <- mask | istate==collapsed.calls$state
	}
	collapsed.calls <- collapsed.calls[mask,]
	if (nrow(collapsed.calls) == 0) {
		stop("No regions to export!")
	}

	## Write to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	if (separate.tracks) {
		## Export peak calls
		bin <- dec2bin(collapsed.calls$state, ndigits=length(multi.hmm$IDs))
		colnames(bin) <- multi.hmm$IDs
		for (icol in 1:ncol(bin)) {
			numsegments <- length(which(bin[,icol]))
			priority <- 52 + 4*icol
			mask <- bin[,icol]
			if (is.null(collapsed.calls$combination)) {
				df <- cbind(collapsed.calls[mask,c('chromosome','start','end','state')], score=rep(0,numsegments), strand=rep(".",numsegments))
			} else {
				df <- cbind(collapsed.calls[mask,c('chromosome','start','end','combination')], score=rep(0,numsegments), strand=rep(".",numsegments))
			}
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$start
			df$thickEnd <- df$end
			RGB <- t(col2rgb(state.colors['modified']))
			RGB <- apply(RGB,1,paste,collapse=",")
			df$itemRgb <- rep(RGB, numsegments)
			if (header) {
				cat(paste0("track name=\"multivariate calls for ",colnames(bin)[icol],"\" description=\"multivariate calls for ",colnames(bin)[icol],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			}
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	} else {
		## Export combinatorial states
		# Generate the colors for each combinatorial state
		colors <- rainbow(numstates)
		RGBs <- t(col2rgb(colors))
		RGBs <- apply(RGBs,1,paste,collapse=",")
		itemRgb <- RGBs[as.integer(factor(as.character(collapsed.calls$state)))]

		# Write to file
		numsegments <- nrow(collapsed.calls)
		if (is.null(collapsed.calls$combination)) {
			df <- cbind(collapsed.calls[,c('chromosome','start','end','state')], score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		} else {
			df <- cbind(collapsed.calls[,c('chromosome','start','end','combination')], score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		}
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		df$thickStart <- df$thickStart - 1
		if (header) {
			cat(paste0("track name=\"combinatorial state\" description=\"multivariate combinatorial states\" visibility=1 itemRgb=On priority=49\n"), file=filename.gz, append=TRUE)
		}
		write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	close(filename.gz)

}


#----------------------------------------------------
# Export read counts from multivariate HMMs
#----------------------------------------------------
exportMultivariateReadCounts <- function(multi.hmm, filename="chromstaR_multivariateReadCounts", header=TRUE) {

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
	filename <- paste0(filename,".wiggle.gz")
	filename.gz <- gzfile(filename, 'w')
	nummod <- length(multi.hmm$IDs)
	readcol <- paste(col2rgb(state.colors['reads']), collapse=',')

	## Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		ID <- multi.hmm$IDs[imod]
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
			write.table(multi.hmm$bins[seqnames(multi.hmm$bins)==chrom]$reads[,imod], file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
}


#=====================================================
# Export binned data
#=====================================================
#' Export genome browser viewable files
#'
#' Export read counts as genome browser viewable file
#'
#' Export read counts from \code{\link{binned.data}} as a file which can be uploaded into a genome browser. Read counts are exported in WIGGLE format (.wiggle.gz).
#'
#' @author Aaron Taudt
#' @param binned.data.list A list of \code{\link{binned.data}} objects or files that contain such objects.
#' @param filename The name of the file that will be written. The ending ".wiggle.gz" for read counts will be appended. Any existing file will be overwritten.
#' @param header A logical indicating whether the output file will have a heading track line (\code{TRUE}) or not (\code{FALSE}).
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportMultivariate}}
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
exportBinnedData <- function(binned.data.list, filename="chromstaR_ReadCounts", header=TRUE) {

	## Load data
	if (is.character(binned.data.list)) {
		message('Loading binned.data from files ...', appendLF=F); ptm <- proc.time()
		binfiles <- binned.data.list
		binned.data.list <- list()
		for (binfile in binfiles) {
			binned.data.list[[binfile]] <- get(load(binfile))
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
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
	filename <- paste0(filename,".wiggle.gz")
	filename.gz <- gzfile(filename, 'w')
	readcol <- paste(col2rgb(state.colors['reads']), collapse=',')

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing binned data ',imod,' / ',nummod,'\r', appendLF=F)
		b <- binned.data.list[[imod]]
		priority <- 50 + 4*imod
		binsize <- width(b[1])
		name <- names(binned.data.list)[imod]
		if (header) {
			cat(paste0('track type=wiggle_0 name="read count for ',name,'" description="read count for ',name,'" visibility=full autoScale=on color=',readcol,' maxHeightPixels=100:50:20 graphType=bar priority=',priority,'\n'), file=filename.gz, append=TRUE)
		}
		# Write read data
		for (chrom in unique(b$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			write.table(mcols(b[b$chromosome==chrom])$reads, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
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
#' @seealso \code{\link{exportUnivariates}}, \code{\link{exportMultivariate}}
#' @export
exportGRanges <- function(gr, trackname, filename="chromstaR_GRanges_regions", header=TRUE) {

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
	filename.gz <- gzfile(filename, 'w')

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	if (header) {
		cat(paste0("track name=\"",trackname,"\" description=\"",trackname,"\" visibility=1 itemRgb=Off\n"), file=filename.gz, append=TRUE)
	}
	if (is.null(gr$score)) {
		gr$score <- 0
	}
	regions <- as.data.frame(gr)[c('chromosome','start','end','score')]
	regions$score <- as.integer(regions$score)
	regions$name <- paste0('region_', 1:nrow(regions))
	regions <- regions[c('chromosome','start','end','name','score')]
	numsegments <- nrow(regions)
	df <- cbind(regions, strand=rep(".",numsegments), thickStart=regions$start, thickEnd=regions$end)
	# Convert from 1-based closed to 0-based half open
	df$start <- df$start - 1
	df$thickStart <- df$thickStart - 1
	if (nrow(df) == 0) {
		warning('No regions in input')
	} else {
		write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

	close(filename.gz)
	message('')

}


