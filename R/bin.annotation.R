bin.annotation = function(annotation.file, reference.genome.file, outputfolder="binned_annotation", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE) {

	# Functions
	get_point_feature = function(coords) {
		out = rep(0,numbins)
		indices = sort(coords %/% binsize + 1)
		rle_indices = rle(indices)
		out[rle_indices$values] = rle_indices$lengths
		return(out)
	}
	get_range_feature = function(starts, ends) {
		out = rep(0,numbins)
		indices = matrix(c(starts%/%binsize+1, ends%/%binsize+1), ncol=2)
		if (length(indices) > 0) {
			for (i1 in 1:nrow(indices)) {
				out[indices[i1,1]:indices[i1,2]] = out[indices[i1,1]:indices[i1,2]] + 1
			}
		}
		return(out)
	}

	# Read in the data
	cat("Reading file",basename(annotation.file),"...")
	annotation = utils::read.table(annotation.file, sep="\t", header=TRUE)
	annotation$txDiff = annotation$txEnd - annotation$txStart
	annotation$cdsDiff = annotation$cdsEnd - annotation$cdsStart
	cat(" done\n")

	# Reference genome
	reference.genome = utils::read.table(reference.genome.file, row.names=1)

	# Do the loop for all binsizes
	for (binsize in binsizes) {
	cat("Binning into binsize",binsize,"\n")

		# Iterate over all chromosomes
		binned.annotation.allchroms = NULL
		chroms.in.annotation = levels(annotation$chrom)
		if (is.null(chromosomes)) {
			chromosomes = chroms.in.annotation
		}
		for (chromosome in chromosomes) {
			# Check if chromosome exists in annotation
			if ( !(chromosome %in% chroms.in.annotation) ) {
				cat(chromosome,"is not in the annotation! Skipped.\n")
				next
			} else if ( !(chromosome %in% rownames(reference.genome)) ) {
				cat(chromosome,"is not in the reference annotation! Skipped.\n")
				next
			}
			cat(chromosome,"              \n")
			## Get annotation for this chromosome
			cat("get annotations...         \r")
			numbins = floor(reference.genome[chromosome,]/binsize)
			iannotation = annotation[annotation$chrom==chromosome, ]
			## Make bin coordinates
			cat("make bin coordinates...    \r")
			bin.chroms = rep(chromosome,numbins)
			bin.starts = seq(from=0, by=binsize, length.out=numbins)
			bin.ends = seq(from=binsize-1, by=binsize, length.out=numbins)
			## Map annotations to their bin
			cat("map annotations to bin...  \r")

			TSS = get_point_feature(iannotation$txStart)
			TES = get_point_feature(iannotation$txEnd)
			TRS = get_range_feature(iannotation$txStart,iannotation$txEnd)

			CSS = get_point_feature(iannotation$cdsStart)
			CES = get_point_feature(iannotation$cdsEnd)
			CRS = get_range_feature(iannotation$cdsStart,iannotation$cdsEnd)


			exonStarts_list = lapply(strsplit(as.character(iannotation$exonStarts), ","), as.integer)
			exonStarts = unlist(exonStarts_list)
			exonEnds_list = lapply(strsplit(as.character(iannotation$exonEnds), ","), as.integer)
			exonEnds = unlist(exonEnds_list)
			ESS = get_point_feature(exonStarts)
			EES = get_point_feature(exonEnds)
			ERS = get_range_feature(exonStarts,exonEnds)

			# Introns
			intronStarts_list = NULL
			intronEnds_list = NULL
			for (i1 in 1:nrow(iannotation)) {
				intron.end = NULL
				intron.start = NULL
				intron.start = iannotation$txStart[i1]
				for (i2 in 1:length(exonStarts_list[[i1]])) {
					intron.end[i2] = exonStarts_list[[i1]][i2] - 1
					intron.start[i2+1] = exonEnds_list[[i1]][i2] + 1
				}
				intron.end[i2+1] = iannotation$txEnd[i1]
				mask = intron.start<=intron.end
				intronStarts_list[[i1]] = intron.start[mask]
				intronEnds_list[[i1]] = intron.end[mask]
			}
			intronStarts = unlist(intronStarts_list)
			intronEnds = unlist(intronEnds_list)
			ISS = get_point_feature(intronStarts)
			IES = get_point_feature(intronEnds)
			IRS = get_range_feature(intronStarts,intronEnds)

			# Protein coding
			mask = iannotation$proteinID!=""
			protein.coding = get_range_feature(iannotation$txStart[mask], iannotation$txEnd[mask])

			# Concatenate
			cat("concatenate...             \r")
			binned.annotation = data.frame(
				chrom = bin.chroms,
				start = bin.starts,
				end = bin.ends,
				TSS = TSS,
				TRS = TRS,
				TES = TES,
				CSS = CSS,
				CRS = CRS,
				CES = CES,
				ESS = ESS,
				ERS = ERS,
				EES = EES,
				ISS = ISS,
				IRS = IRS,
				IES = IES,
				PC = protein.coding
			)
			
			if (separate.chroms==TRUE) {
				if (save.as.RData==TRUE) {
					# Print to file
					filename = paste(basename(annotation.file),"_binsize_",binsize,"_",chromosome,".RData", sep="")
					cat("save...                    \r")
					save(binned.annotation, file=file.path(outputfolder,filename) )
				} else {
					cat("                          \r")
					return(binned.annotation)
				}
			} else {
				binned.annotation.allchroms[[length(binned.annotation.allchroms)+1]] = binned.annotation
			}
			cat("                          \r")

		}
		if (separate.chroms!=TRUE) {
			cat("Concatenating chromosomes ...")
			binned.annotation.allchroms = do.call("rbind",binned.annotation.allchroms)
			cat(" done\n")
			if (save.as.RData==TRUE) {
				# Print to file
				filename = paste(basename(annotation.file),"_binsize_",binsize,".RData", sep="")
				cat("Saving to file ...")
				save(binned.annotation.allchroms, file=file.path(outputfolder,filename) )
				cat(" done\n")
			} else {
				return(binned.annotation.allchroms)
			}
		}

	}

}


