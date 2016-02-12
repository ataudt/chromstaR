#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
#'
#' Convert aligned reads or signal values from various file formats into read counts in equidistant windows (bins). bam2binned() and bed2binned() use 'GenomicRanges::countOverlaps()' to calculate the read counts. Therefore, with small binsizes and large read lengths, one read fragment will often overlap more than one bin.
#' bedGraph2binned() sets the maximum signal value in a bin as value for that bin.
#'
#' @author Aaron Taudt
#' @examples
#'## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue
#'bedfile <- system.file(
#'      "extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz",
#'      package="chromstaR")
#'## Bin the BED files into bin size 1000bp
#'binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000,
#'                          save.as.RData=FALSE)
#'## Get the number of bins
#'length(binned.data)
#' @seealso \link{binned.data}
#' @name binning
NULL

#' @describeIn binning Bin reads in BAM format.
#' @param bamfile A file in BAM format.
#' @inheritParams align2binned
#' @param bamindex BAM index file. Can be specified without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @export
bam2binned <- function(bamfile, bamindex=bamfile, pairedEndReads=FALSE, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=FALSE, min.mapq=10, downsample.to.reads=NULL, remove.duplicate.reads=FALSE, max.fragment.width=1000) {
	return(align2binned(bamfile, format="bam", bamindex=bamindex, pairedEndReads=pairedEndReads, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData, min.mapq=min.mapq, downsample.to.reads=downsample.to.reads, remove.duplicate.reads=remove.duplicate.reads, max.fragment.width=max.fragment.width))
}

#' @describeIn binning Bin reads in BED format.
#' @param bedfile A file in BED format.
#' @inheritParams align2binned
#' @export
bed2binned <- function(bedfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=FALSE, downsample.to.reads=NULL) {
	return(align2binned(bedfile, format="bed", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData, downsample.to.reads=downsample.to.reads))
}

#' @describeIn binning Bin reads in bedGraph format.
#' @param bedGraphfile A file in bedGraph format.
#' @inheritParams align2binned
#' @export
bedGraph2binned <- function(bedGraphfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=FALSE, downsample.to.reads=NULL) {
	return(align2binned(bedGraphfile, format="bedGraph", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData, downsample.to.reads=downsample.to.reads))
}

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
#'
#' @param file A file with aligned reads.
#' @param format One of \code{c('bam', 'bed', 'bedGraph')}.
#' @param assembly An assembly to specify the chromosome lengths. One of \code{c('hg19','hg18')}.
#' @param bamindex Index file if \code{format='bam'} with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param chrom.length.file A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.
#' @param outputfolder Folder to which the binned data will be saved if \code{save.as.RData=TRUE}. If the specified folder does not exist, it will be created.
#' @param binsizes A vector with integer values which will be used for the binning. If more than one value is given, output files will be produced for each bin size.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param save.as.RData If set to \code{TRUE}, output will be written to file. The filename is generated from the input file and bin sizes. If set to \code{FALSE}, no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param downsample.to.reads Downsample the input data to the specified number of reads.
#' @param remove.duplicate.reads Logical indicating whether duplicate reads should be removed. NOTE: Duplicate removal only works if your duplicate reads are properly flagged in the BAM file.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @return The function produces a \code{\link{GRanges}} object with one meta data column 'reads' that contains the read count. This binned data will be either written to file (\code{save.as.RData=FALSE}) or given as return value (\code{save.as.RData=FALSE}).
align2binned <- function(file, format, assembly, bamindex=file, pairedEndReads=FALSE, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=FALSE, min.mapq=10, downsample.to.reads=NULL, remove.duplicate.reads=FALSE, max.fragment.width=1000) {

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}

	### Read in the data
	message("Reading file ",basename(file)," ...", appendLF=F); ptm <- proc.time()
	## BED (0-based)
	if (format == "bed") {
		if (!is.null(chrom.length.file)) {
			# File with chromosome lengths (1-based)
			chrom.lengths.df <- read.table(chrom.length.file)
			chrom.lengths <- chrom.lengths.df[,2]
			names(chrom.lengths) <- chrom.lengths.df[,1]
		} else {
			chrom.lengths <- getChromLengths(assembly=assembly)
		}
		# File with reads, determine classes first for faster import (0-based)
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:3] <- classes.in.bed[1:3]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)))	# start+1 to go from [0,x) -> [1,x]
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	## BAM (1-based)
	} else if (format == "bam") {
		data <- suppressMessages( bam2GRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, keep.duplicate.reads=!remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width) )
	## BEDGraph (0-based)
	} else if (format == "bedGraph") {
		if (!is.null(chrom.length.file)) {
			# File with chromosome lengths (1-based)
			chrom.lengths.df <- read.table(chrom.length.file)
			chrom.lengths <- chrom.lengths.df[,2]
			names(chrom.lengths) <- chrom.lengths.df[,1]
		} else {
			chrom.lengths <- getChromLengths(assembly=assembly)
		}
		# File with reads, determine classes first for faster import
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:4] <- classes.in.bed[1:4]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)), signal=data[,4])	# start+1 to go from [0,x) -> [1,x]
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Select chromosomes to bin
	if (is.null(chromosomes)) {
		chromosomes <- seqlevels(data)
	}
	chroms2use <- intersect(chromosomes, seqlevels(data))
 
	### Downsampling
	if (!is.null(downsample.to.reads)) {
		if (downsample.to.reads <= length(data)) {
			message("downsampling to ", downsample.to.reads, " reads ...", appendLF=F); ptm <- proc.time()
			data <- data[sort(sample(1:length(data), size=downsample.to.reads, replace=FALSE))]
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
			warning("No downsampling done because there are only ", length(data), " reads in the data.")
		}
	}

	### Do the loop for all binsizes
	for (binsize in binsizes) {
		message("Binning into bin size ",binsize)

		### Iterate over all chromosomes
		message("  binning genome ...", appendLF=F); ptm <- proc.time()
		binned.data <- GenomicRanges::GRangesList()
		for (chromosome in chroms2use) {
			## Check last incomplete bin
			incomplete.bin <- seqlengths(data)[chromosome] %% binsize > 0
			if (incomplete.bin) {
				numbins <- floor(seqlengths(data)[chromosome]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
			} else {
				numbins <- seqlengths(data)[chromosome]/binsize
			}
			if (numbins == 0) {
				warning("chromosome ",chromosome," is smaller than binsize. Skipped.")
				next
			}
			## Initialize vectors
			chroms <- rep(chromosome,numbins)
			reads <- rep(0,numbins)
			start <- seq(from=1, by=binsize, length.out=numbins)
			end <- seq(from=binsize, by=binsize, length.out=numbins)
# 			end[length(end)] <- seqlengths(data)[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired

			## Create binned chromosome as GRanges object
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbins),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbins)
							)
			seqlengths(i.binned.data) <- seqlengths(data)[chromosome]

			suppressWarnings(
				binned.data[[chromosome]] <- i.binned.data
			)
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		binned.data <- unlist(binned.data)
		names(binned.data) <- NULL

		## Count overlaps
		message("  counting overlaps ...", appendLF=F); ptm <- proc.time()
		if (format=="bam" | format=="bed") {
			reads <- GenomicRanges::countOverlaps(binned.data, data)
		} else if (format=="bedGraph") {
			# Take the max value from all regions that fall into / overlap a given bin as read count
			midx <- as.matrix(findOverlaps(binned.data, data))
			reads <- rep(0,length(binned.data))
			signal <- mcols(data)$signal
			rle <- rle(midx[,1])
			read.idx <- rle$values
			max.idx <- cumsum(rle$lengths)
			maxvalues <- rep(NA, length(read.idx))
			maxvalues[1] <- max(signal[midx[1:(max.idx[1]),2]])
			for (i1 in 2:length(read.idx)) {
				maxvalues[i1] <- max(signal[midx[(max.idx[i1-1]+1):(max.idx[i1]),2]])
			}
			reads[read.idx] <- maxvalues
		}
		mcols(binned.data)$reads <- reads
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		if (length(binned.data) == 0) {
			warning(paste0("The bin size of ",binsize," is larger than any of the chromosomes."))
			return(NULL)
		}

		if (save.as.RData==TRUE) {
			# Print to file
			filename <- paste0(basename(file),"_binsize",format(binsize, scientific=F, trim=T),".RData")
			message("Saving to file ...", appendLF=F); ptm <- proc.time()
			save(binned.data, file=file.path(outputfolder,filename) )
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
			return(binned.data)
		}

	}

}


#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param file A file in BAM format.
#' @param bamindex BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be imported, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param keep.duplicate.reads A logical indicating whether or not duplicate reads should be kept.
#' @param min.mapq Minimum mapping quality when importing from BAM files. Set \code{min.mapq=NULL} to keep all reads.
#' @param max.fragment.width Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairsFromBam readGAlignmentsFromBam first
bam2GRanges <- function(file, bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, keep.duplicate.reads=TRUE, min.mapq=10, max.fragment.width=1000) {

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		bamindex.own <- Rsamtools::indexBam(file)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
	}
	file.header <- Rsamtools::scanBamHeader(file)[[1]]
	chrom.lengths <- file.header$targets
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data.')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}
	## Import the file into GRanges
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (keep.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		}
	}
	## Filter by mapping quality
	if (pairedEndReads) {
		message("Converting to GRanges ...", appendLF=FALSE); ptm <- proc.time()
		data <- as(data.raw, 'GRanges') # treat as one fragment
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		message("Filtering reads ...", appendLF=FALSE); ptm <- proc.time()
		if (!is.null(min.mapq)) {
			mapq.first <- mcols(GenomicAlignments::first(data.raw))$mapq
			mapq.last <- mcols(GenomicAlignments::last(data.raw))$mapq
			mapq.mask <- mapq.first >= min.mapq & mapq.last >= min.mapq
			if (any(is.na(mapq.mask))) {
				warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
			}
			data <- data[which(mapq.mask)]
			# Filter out too long fragments
			data <- data[width(data)<=max.fragment.width]
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	} else {
		message("Converting to GRanges ...", appendLF=FALSE); ptm <- proc.time()
		data <- as(data.raw, 'GRanges')
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		message("Filtering reads ...", appendLF=FALSE); ptm <- proc.time()
		if (!is.null(min.mapq)) {
			if (any(is.na(mcols(data)$mapq))) {
				warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
				mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
			}
			data <- data[mcols(data)$mapq >= min.mapq]
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	}
	return(data)

}

