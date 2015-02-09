# ==============================================
# Write color-coded tracks with states from HMMs
# ==============================================
export.unihmm2bed <- function(hmm.list, only.modified=TRUE, filename="view_me_in_genome_browser") {

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
		priority <- 51 + 3*imod
		cat(paste0("track name=\"univariate calls for ",hmm$ID,"\" description=\"univariate calls for ",hmm$ID,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
		collapsed.calls <- as.data.frame(hmm.gr)[c('chromosome','start','end','state')]
		itemRgb <- RGBs[as.character(collapsed.calls$state)]
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		df$thickStart <- df$thickStart - 1
		if (only.modified) {
			df <- df[df$state=='modified',]
		}
		if (nrow(df) == 0) {
			warning('hmm ',imod,' does not contain any \'modified\' calls')
		} else {
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}
	close(filename.gz)
	message('')

}


# =============================
# Write signal tracks from HMMs
# =============================
export.unihmm2wiggle <- function(hmm.list, filename="view_me_in_genome_browser") {

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

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		hmm <- hmm.list[[imod]]
		hmm.gr <- hmm.grl[[imod]]
		priority <- 50 + 3*imod
		binsize <- width(hmm.gr[1])
		cat(paste0("track type=wiggle_0 name=\"read count for ",hmm$ID,"\" description=\"read count for ",hmm$ID,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
		# Write read data
		for (chrom in unique(hmm.gr$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			write.table(mcols(hmm.gr[hmm.gr$chromosome==chrom])$reads, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
}


# ===============================================================
# Write color-coded tracks with multivariate combinatorial states
# ===============================================================
export.multihmm2bed <- function(multi.hmm, separate.tracks=TRUE, exclude.states=0, include.states=NULL, filename="view_me_in_genome_browser") {

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
	collapsed.calls <- as.data.frame(insertchr(multi.hmm$segments))[,c(7,2:3,6)]
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
		bin <- dec2bin(collapsed.calls$state, ndigits=length(multi.hmm$IDs))
		colnames(bin) <- multi.hmm$IDs
		for (icol in 1:ncol(bin)) {
			numsegments <- length(which(bin[,icol]))
			priority <- 52 + 3*icol
			df <- cbind(collapsed.calls[bin[,icol],], score=rep(0,numsegments), strand=rep(".",numsegments))
			# Convert from 1-based closed to 0-based half open
			df$start <- df$start - 1
			df$thickStart <- df$start
			df$thickEnd <- df$end
			RGB <- t(col2rgb(state.colors['modified']))
			RGB <- apply(RGB,1,paste,collapse=",")
			df$itemRgb <- rep(RGB, numsegments)
			cat(paste0("track name=\"multivariate calls for ",colnames(bin)[icol],"\" description=\"multivariate calls for ",colnames(bin)[icol],"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	} else {
		# Generate the colors for each combinatorial state
		colors <- colors()[grep(colors(), pattern="white|grey|gray|snow|aliceblue|azure", invert=T)]
		step <- length(colors) %/% numstates
		colors <- colors[seq(1,by=step,length=numstates)]
		RGBs <- t(col2rgb(colors))
		RGBs <- apply(RGBs,1,paste,collapse=",")
		itemRgb <- RGBs[as.integer(factor(as.character(collapsed.calls$state)))]

		# Write to file
		numsegments <- nrow(collapsed.calls)
		df <- cbind(collapsed.calls, score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=collapsed.calls$start, thickEnd=collapsed.calls$end, itemRgb=itemRgb)
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		df$thickStart <- df$thickStart - 1
		cat(paste0("track name=\"combinatorial state\" description=\"multivariate combinatorial states\" visibility=1 itemRgb=On priority=49\n"), file=filename.gz, append=TRUE)
		write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	close(filename.gz)

}


# ===================================================
# Write tracks with read counts from multivariate HMM
# ===================================================
export.multihmm2wiggle <- function(multi.hmm, filename="view_me_in_genome_browser") {

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

	## Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing hmm ',imod,' / ',nummod,'\r', appendLF=F)
		ID <- multi.hmm$IDs[imod]
		priority <- 50 + 3*imod
		binsize <- width(multi.hmm$bins[1])
		cat(paste0("track type=wiggle_0 name=\"read count for ",ID,"\" description=\"read count for ",ID,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
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


# =====================
# Export GRanges to bed
# =====================
export.GRanges2bed <- function(gr, filename="view_me_in_genome_browser", separate.tracks=FALSE) {

	# Variables
	filename <- paste0(filename,".bed.gz")
	filename.gz <- gzfile(filename, 'w')
	numstates <- length(unique(gr$state))

	# Generate the colors for each combinatorial state
	colors <- colors()[grep(colors(), pattern="white|grey|gray|snow|aliceblue|azure", invert=T)]
	step <- length(colors) %/% numstates
	colors <- colors[seq(1,by=step,length=numstates)]
	RGBs <- t(col2rgb(colors))
	RGBs <- apply(RGBs,1,paste,collapse=",")
	itemRgb <- RGBs[as.integer(factor(gr$state))]

	# Write to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	numsegments <- length(gr)
	df <- cbind(as.data.frame(gr)[,c(1:3,6)], score=rep(0,numsegments), strand=rep(".",numsegments), thickStart=start(ranges(gr)), thickEnd=end(ranges(gr)), itemRgb=itemRgb)
	# Adjust coordinates to BED format
	df$start <- df$start - 1

	if (separate.tracks) {
		for (istate in unique(gr$state)) {
			priority <- 100 + which(istate==unique(gr$state))
			cat(paste0("track name=\"combinatorial state ",istate,"\" description=\"combinatorial state ",istate,"\" visibility=1 itemRgb=On priority=",priority,"\n"), file=filename.gz, append=TRUE)
			write.table(format(df[df$name==istate,], scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	} else {
		cat(paste0("track name=\"combinatorial states\" description=\"combinatorial states\" visibility=1 itemRgb=On priority=100\n"), file=filename.gz, append=TRUE)
		write.table(format(df, scientific=FALSE), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	close(filename.gz)

}


# =================================
# Export binned GRanges to bedGraph
# =================================
export.binned2bedGraph <- function(binned.data.list, filename="view_me_in_genome_browser") {

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
	nummod <- length(binned.data.list)
	filename <- paste0(filename,".bedGraph.gz")
	filename.gz <- gzfile(filename, 'w')

	## Load data
	for (imod in 1:nummod) {
		if (class(binned.data.list[[imod]])!='GRanges') {
			binned.data.list[[imod]] <- get(load(binned.data.list[[imod]]))
			if (class(binned.data.list[[imod]])!='GRanges') stop("argument 'binned.data.list' expects a list with GRanges objects with meta column 'reads'")
		}
	}

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing binned.data ',imod,' / ',nummod,'\r', appendLF=F)
		gr <- binned.data.list[[imod]]
		df <- as.data.frame(gr)[,c(1,2,3,6)]
		df <- df[df$reads>0,]
		# Convert from 1-based closed to 0-based half open
		df$start <- df$start - 1
		priority <- 50 + 3*imod
		name <- names(binned.data.list)[imod]
		cat(paste0("track type=bedGraph name=\"read count for ",name,"\" description=\"read count for ",name,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
		# Write read data
		write.table(format(df, scientific=F), file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE, quote=F)
	}
	close(filename.gz)
	message('')

}


# ====================================
# Write signal tracks from binned data
# ====================================
export.binned2wiggle <- function(binned.data.list, filename="view_me_in_genome_browser") {

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

	# Write first line to file
	message('writing to file ',filename)
	cat("", file=filename.gz)
	
	### Write every model to file ###
	for (imod in 1:nummod) {
		message('writing binned data ',imod,' / ',nummod,'\r', appendLF=F)
		b <- binned.data.list[[imod]]
		priority <- 50 + 3*imod
		binsize <- width(b[1])
		name <- names(binned.data.list)[imod]
		cat(paste0("track type=wiggle_0 name=\"read count for ",name,"\" description=\"read count for ",name,"\" visibility=full autoScale=on color=90,90,90 maxHeightPixels=100:50:20 graphType=bar priority=",priority,"\n"), file=filename.gz, append=TRUE)
		# Write read data
		for (chrom in unique(b$chromosome)) {
			cat(paste0("fixedStep chrom=",chrom," start=1 step=",binsize," span=",binsize,"\n"), file=filename.gz, append=TRUE)
			write.table(mcols(b[b$chromosome==chrom])$reads, file=filename.gz, append=TRUE, row.names=FALSE, col.names=FALSE)
		}
	}
	close(filename.gz)
	message('')
}


