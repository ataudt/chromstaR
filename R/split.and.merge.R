split.per.chrom <- function(uni.hmm, filename=NULL) {

	## Check user input
	if (check.univariate.model(uni.hmm)!=0) {
		uni.hmm <- get(load(uni.hmm))
		if (check.univariate.model(uni.hmm)!=0) stop("argument 'uni.hmm' expects a univariate HMM or a file that contains a univariate HMM")
	}

	## Construct new HMM per chromosome
	split.hmm <- uni.hmm
	split.hmm[c('bins','segments')] <- NULL
	split.bins <- split(uni.hmm$bins, seqnames(uni.hmm$bins))
	split.segments <- split(uni.hmm$segments, seqnames(uni.hmm$segments))

	hmm.list <- list()
	for (chrom in seqlevels(uni.hmm$bins)) {
		chr.hmm <- uni.hmm
		chr.hmm$bins <- split.bins[[chrom]]
		chr.hmm$segments <- split.segments[[chrom]]
		hmm.list[[chrom]] <- chr.hmm
		if (!is.null(filename)) {
			save(chr.hmm, file=paste0(filename,'_chr_',chrom,'.RData'))
		}
	}
	if (is.null(filename)) {
		return(hmm.list)
	}

}

merge.chroms <- function(multi.hmm.list, filename=NULL) {

	## Check user input
	multi.hmm.list <- chromstar:::loadMultiHmmsFromFiles(multi.hmm.list)
		
	## Check if posteriors are present everywhere and adjust number of columns if necessary
	num.models <- length(multi.hmm.list)
	num.tracks.per.model <- unlist(lapply(multi.hmm.list, function(x) { length(levels(x$bins$state)) }))
	num.tracks <- max(num.tracks.per.model)
	post.colnames <- paste0("P(",levels(multi.hmm.list[[which.max(num.tracks.per.model)]]$bins$state),")")
	post.present <- Reduce('|', unlist(lapply(multi.hmm.list, function(x) { !is.null(x$bins$posteriors) })))
	
	## Construct new HMM
	cat("concatenating HMMs\n")
	bins <- list()	# do not use GRangesList() because it copies the whole list each time an element is added
	segments <- list()
	for (i1 in 1:num.models) {
		cat(" ",i1," of ",num.models,"        \r")
		hmm <- multi.hmm.list[[1]]	# select always first because we remove it at the end of the loop
		if (post.present) {
			if (length(levels(hmm)) < num.tracks) {
				missing.post.colnames <- setdiff(post.colnames, colnames(hmm$bins$posteriors))
				missing.posteriors <- matrix(0, ncol=length(missing.post.colnames), nrow=nrow(hmm$bins$posteriors))
				colnames(missing.posteriors) <- missing.post.colnames
				posteriors <- cbind(hmm$bins$posteriors, missing.posteriors)
				remove(missing.posteriors)
				hmm$bins$posteriors <- posteriors[,post.colnames]
				remove(posteriors)
			}
		}
		bins[[i1]] <- hmm$bins
		segments[[i1]] <- hmm$segments
		# Remove current HMM to save memory
		if (i1 < num.models) remove(hmm)	# remove it because otherwise R will make a copy when we NULL the underlying reference (multi.hmm.list[[1]])
		multi.hmm.list[[1]] <- NULL
	}
	cat("\n")
# 	if (!post.present) {
		cat("merging ..."); ptm <- proc.time()
		bins <- do.call('c', bins)	# this can be too memory intensive if posteriors are present
		segments <- do.call('c', segments)
		time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))
# 	} else {
# 		cat("temporarily storing posteriors ..."); ptm <- proc.time()
# 		rm(.Random.seed, envir=globalenv())
# 		randomString <- paste(sample(c(0:9, letters, LETTERS), 15, replace=T), collapse="")
# 		tempfile <- paste0('temp_R_chromstaR_merge.chroms_posteriors_',randomString)
# # 		tempfile.gz <- gzfile(tempfile, 'w')
# 		tempfile.gz <- tempfile
# 		for (i1 in 1:length(bins)) {
# 			if (i1 == 1) {
# 				write.table(bins[[i1]]$posteriors, row.names=F, col.names=F, file=tempfile.gz, append=F)
# 			} else {
# 				write.table(bins[[i1]]$posteriors, row.names=F, col.names=F, file=tempfile.gz, append=T)
# 			}
# 		}
# # 		close(tempfile.gz)
# 		time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))
# 		cat("merging ..."); ptm <- proc.time()
# 		bins <- do.call('c', bins)
# 		segments <- do.call('c', segments)
# 		time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))
# 		cat("loading posteriors ..."); ptm <- proc.time()
# 		tab5rows <- read.table(tempfile, nrows=5)
# 		classes.in.file <- sapply(tab5rows, class)
# 		bins$posteriors <- as.matrix(read.table(tempfile, header=F, colClasses=classes.in.file))
# 		colnames(bins$posteriors) <- post.colnames
# 		file.remove(tempfile)
# 		time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))
# 	}
	multi.hmm <- hmm
	multi.hmm$bins <- bins
	multi.hmm$segments <- segments
	cat("calculating weights ..."); ptm <- proc.time()
	if (post.present) {
		multi.hmm$weights <- apply(multi.hmm$bins$posteriors, 2, mean)
		names(multi.hmm$weights) <- levels(multi.hmm$bins$state)
	} else {
		multi.hmm$weights <- table(multi.hmm$bins$state) / length(multi.hmm$bins)
	}
	time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))

	if (is.null(filename)) {
		return(multi.hmm)
	} else {
		cat("writing to file",filename,"..."); ptm <- proc.time()
		save(multi.hmm, file=filename)
		time <- proc.time() - ptm; cat(paste0(" ",round(time[3],2),"s\n"))
	}

}
