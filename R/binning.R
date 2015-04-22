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
#' @param bamindex BAM index file. Can be specified without the .bai ending.
#' @export
bam2binned <- function(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=TRUE) {
	return(align2binned(bamfile, format="bam", index=bamindex, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData))
}

#' @describeIn binning Bin reads in BED format.
#' @param bedfile A file in BED format.
#' @inheritParams align2binned
#' @export
bed2binned <- function(bedfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=TRUE) {
	return(align2binned(bedfile, format="bed", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData))
}

#' @describeIn binning Bin reads in bedGraph format.
#' @param bedGraphfile A file in bedGraph format.
#' @inheritParams align2binned
#' @export
bedGraph2binned <- function(bedGraphfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=TRUE) {
	return(align2binned(bedGraphfile, format="bedGraph", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, save.as.RData=save.as.RData))
}

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
#'
#' @param file A file with aligned reads.
#' @param format One of \code{c('bam', 'bed', 'bedGraph')}.
#' @param assembly An assembly to specify the chromosome lengths. One of \code{c('hg19','hg18')}.
#' @param index Index file if \code{format='bam'} with or without the .bai ending.
#' @param chrom.length.file A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.
#' @param outputfolder Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.
#' @param binsizes A vector with integer values which will be used for the binning. If more than one value is given, output files will be produced for each bin size.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param save.as.RData If set to \code{FALSE}, no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case.
#' @return The function produces a \link{GRanges} object with one meta data column 'reads' that contains the read count. This binned data will be either written to file (\code{save.as.RData=TRUE}) or given as return value (\code{save.as.RData=FALSE}).
#' @import Rsamtools
#' @import GenomicAlignments
#' @import BSgenome
align2binned <- function(file, format, assembly, index=file, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, save.as.RData=TRUE) {

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
		chroms.in.data <- seqlevels(data)
	## BAM (1-based)
	} else if (format == "bam") {
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
		gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
		data <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(which=range(gr)))
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
		chroms.in.data <- seqlevels(data)
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	## Select chromosomes to bin
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}
	diff <- setdiff(chromosomes, names(chrom.lengths))
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because no lengths could be found.'))
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	chroms2use <- intersect(chroms2use, names(chrom.lengths))
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data.')
	}
 
	### Do the loop for all binsizes
	for (binsize in binsizes) {
		message("Binning into bin size ",binsize)

		### Iterate over all chromosomes
		message("  binning genome ...", appendLF=F); ptm <- proc.time()
		binned.data <- GenomicRanges::GRangesList()
		for (chromosome in chroms2use) {
			## Check last incomplete bin
			incomplete.bin <- chrom.lengths[chromosome] %% binsize > 0
			if (incomplete.bin) {
				numbins <- floor(chrom.lengths[chromosome]/binsize)	# floor: we don't want incomplete bins, ceiling: we want incomplete bins at the end
			} else {
				numbins <- chrom.lengths[chromosome]/binsize
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
# 			end[length(end)] <- chrom.lengths[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired

			## Create binned chromosome as GRanges object
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbins),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbins)
							)
			seqlengths(i.binned.data) <- chrom.lengths[chromosome]

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
			filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F, trim=T),".RData")
			message("Saving to file ...", appendLF=F); ptm <- proc.time()
			save(binned.data, file=file.path(outputfolder,filename) )
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
			return(binned.data)
		}

	}

}
