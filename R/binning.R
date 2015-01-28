bedGraph2binned <- function(bedGraphfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bedGraphfile, format="bedGraph", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bam2binned <- function(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bamfile, format="bam", index=bamindex, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bed2binned <- function(bedfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE) {
	return(align2binned(bedfile, format="bed", assembly=assembly, chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

align2binned <- function(file, format, assembly, index=file, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE) {

	## Check user input
	if (save.as.RData==FALSE) {
		separate.chroms=FALSE
	}

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}

	### Read in the data
	## BED (0-based)
	if (format == "bed") {
		message("Reading file ",basename(file)," ...", appendLF=F)
		if (!is.null(chrom.length.file)) {
			# File with chromosome lengths (1-based)
			chrom.lengths.df <- read.table(chrom.length.file)
			chrom.lengths <- chrom.lengths.df[,2]
			names(chrom.lengths) <- chrom.lengths.df[,1]
		} else {
			chrom.lengths <- get.chrom.lengths(assembly=assembly)
		}
		# File with reads, determine classes first for faster import (0-based)
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:3] <- classes.in.bed[1:3]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GenomicRanges::GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2]+1, end=data[,3]), strand=Rle(strand("*"), nrow(data)))	# start+1 to go from [0,x) -> [1,x]
		tC <- tryCatch({
			seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		}, warning = function(war) {
			suppressWarnings(seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))]))
			offending.chroms <- unique(names(which(end(data) > seqlengths(data)[as.character(seqnames(data))])))
			if (war$message=="'ranges' contains values outside of sequence bounds. See ?trim to subset ranges.") {
				warning(paste0("File \"",file,"\" contains reads outside of sequence bounds. The offending chromosomes were \"",paste(offending.chroms, collapse=', '),"\". Please consider using a different reference assembly."))
			} else {
				print(war)
			}
		})
		chroms.in.data <- seqlevels(data)
	## BAM (1-based)
	} else if (format == "bam") {
		message("Reading header of ",basename(file)," ...", appendLF=F)
		file.header <- Rsamtools::scanBamHeader(file)[[1]]
		chrom.lengths <- file.header$targets
		chroms.in.data <- names(chrom.lengths)
	## BEDGraph (0-based)
	} else if (format == "bedGraph") {
		message("Reading file ",basename(file)," ...", appendLF=F)
		if (!is.null(chrom.length.file)) {
			# File with chromosome lengths (1-based)
			chrom.lengths.df <- read.table(chrom.length.file)
			chrom.lengths <- chrom.lengths.df[,2]
			names(chrom.lengths) <- chrom.lengths.df[,1]
		} else {
			chrom.lengths <- get.chrom.lengths(assembly=assembly)
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
	message(" done")

 
	### Do the loop for all binsizes
	for (binsize in binsizes) {
		message("Binning into binsize ",binsize)

		### Iterate over all chromosomes
		binned.data <- GenomicRanges::GRanges()
		if (is.null(chromosomes)) {
			chromosomes <- chroms.in.data
		}
		for (chromosome in chromosomes) {
			## Check if chromosome exists in data
			if ( !(chromosome %in% chroms.in.data) ) {
				warning("Skipped chromosome ",chromosome,", not in the data!")
				next
			} else if ( !(chromosome %in% names(chrom.lengths)) ) {
				warning("Skipped chromosome ",chromosome,", no length found!")
				next
			}
			message(chromosome,"                              ")
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
			message("initialize vectors...\r", appendLF=F)
			chroms <- rep(chromosome,numbins)
			reads <- rep(0,numbins)
			start <- seq(from=1, by=binsize, length.out=numbins)
			end <- seq(from=binsize, by=binsize, length.out=numbins)
# 			end[length(end)] <- chrom.lengths[chromosome] # last ending coordinate is size of chromosome, only if incomplete bins are desired

			## Create binned chromosome as GRanges object
			message("creating GRanges container...            \r", appendLF=F)
			i.binned.data <- GenomicRanges::GRanges(seqnames = Rle(chromosome, numbins),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbins)
							)
			seqlengths(i.binned.data) <- chrom.lengths[chromosome]

			if (format=="bam") {
				message("reading reads from file...               \r", appendLF=F)
				data <- GenomicAlignments::readGAlignmentsFromBam(file, index=index, param=ScanBamParam(what=c("pos"),which=range(i.binned.data),flag=scanBamFlag(isDuplicate=F)))
			}

			## Count overlaps
			message("counting overlaps...                     \r", appendLF=F)
			if (format=="bam" | format=="bed") {
				reads <- GenomicRanges::countOverlaps(i.binned.data, data[seqnames(data)==chromosome])
			} else if (format=="bedGraph") {
				# Take the max value from all regions that fall into / overlap a given bin as read count
				midx <- as.matrix(findOverlaps(i.binned.data, data[seqnames(data)==chromosome]))
				reads <- rep(0,length(i.binned.data))
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
			
			## Concatenate
			message("concatenate...                           \r", appendLF=F)
			mcols(i.binned.data)$reads <- reads

			if (separate.chroms==TRUE) {
				binned.data <- i.binned.data
				if (save.as.RData==TRUE) {
					## Print to file
					filename <- paste(basename(file),"_binsize_",format(binsize, scientific=F),"_",chromosome,".RData", sep="")
					message("save...                                  \r", appendLF=F)
					save(binned.data, file=file.path(outputfolder,filename) )
				} else {
					message("                                         \r", appendLF=F)
					return(binned.data)
				}
			} else {
				binned.data <- suppressWarnings(BiocGenerics::append(binned.data, i.binned.data))
			}
			message("                                         \r", appendLF=F)

		}
		if (separate.chroms==FALSE) {
			if (save.as.RData==TRUE) {
				# Print to file
				filename <- paste0(basename(file),"_binsize_",format(binsize, scientific=F),".RData")
				message("Saving to file ...", appendLF=F)
				save(binned.data, file=file.path(outputfolder,filename) )
				message(" done")
			} else {
				return(binned.data)
			}
		}

	}

}
