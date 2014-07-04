simulate.univariate = function(coordinates, transition, emission, initial=1) {

	# Calculate some variables
	numstates = ncol(transition)
	if (numstates!=3) {
		stop("The transition matrix is expected to have 3 columns and 3 rows")
	}
	numbins = nrow(coordinates)

	## Make state vector from transition matrix
	# Generate sample of random numbers
	rsample = runif(numbins,0,1)
	# Integrate transition matrix by row and add -1 in front
	cumtransition = cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
	# Generate the state vector by going through each state
	cat("Generating states from transition matrix...")
	states = matrix(rep(NA,numstates*numbins), ncol=numstates)
	for (irow in 1:numstates) {
		states[,irow] = findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
	}
	statevec = rep(NA,numbins)
	statevec[1] = initial
	for (ibin in 2:numbins) {
		statevec[ibin] = states[ibin,statevec[ibin-1]]
	}
	cat(" done\n")

	## Make reads from state vector and distributions
	# Generate the read counts by drawing from distribution
	cat("Generating reads from emission parameters and states...")
	reads = rep(NA, numbins)
	numbins.in.state = aggregate(rep(1,length(statevec)), list(state=statevec), sum)
	reads[statevec==1] = 0
	if (!is.na(numbins.in.state[2,'x'])) {
		reads[statevec==2] = rnbinom(numbins.in.state[2,'x'], size=emission[2,'size'], prob=emission[2,'prob'])
	}
	if (!is.na(numbins.in.state[3,'x'])) {
		reads[statevec==3] = rnbinom(numbins.in.state[3,'x'], size=emission[3,'size'], prob=emission[3,'prob'])
	}
	cat(" done\n")

	# Return the output
	out = list(coordinates = coordinates,
				states = statevec,
				reads = reads,
				transition = transition,
				emission = emission
				)
	return(out)

}



simulate.multivariate = function(coordinates, transition, emissions, weights, sigma, use.states, initial=1) {

	lib = require(mvtnorm)
	if (lib == FALSE) {
		install.packages("mvtnorm")
		library(mvtnorm)
	}

	# Calculate some variables
	numstates = ncol(transition)
	numbins = nrow(coordinates)
	nummod = length(emissions)

	## Make state vector from transition matrix
	# Generate sample of random numbers
	rsample = runif(numbins,0,1)
	# Integrate transition matrix by row and add -1 in front
	cumtransition = cbind(rep(-1,numstates), t(apply(transition, 1, cumsum)))
	# Generate the state vector by going through each state
	cat("Generating states from transition matrix...")
	states = matrix(rep(NA,numstates*numbins), ncol=numstates)
	for (irow in 1:numstates) {
		states[,irow] = findInterval(rsample, cumtransition[irow,], rightmost.closed=TRUE)
	}
	statevec = rep(NA,numbins)
	statevec[1] = initial
	for (ibin in 2:numbins) {
		statevec[ibin] = states[ibin,statevec[ibin-1]]
	}
	# Replace the states by combinatorial states
	statevec = use.states[statevec]
	cat(" done\n")

	## Make reads from state vector and emission distributions
	reads = matrix(rep(NA, nummod*numbins), ncol=nummod)

	sizes = unlist(lapply(emissions, "[", 2:3, 'size'))
	probs = unlist(lapply(emissions, "[", 2:3, 'prob'))
	ws = unlist(lapply(weights, "[", 1))
	for (istate in use.states) {
		cat("Generating reads for state ",istate,"\n")
		i = which(use.states==istate)
		
		# Convert istate to binary representation
		binary_state = rev(as.integer(intToBits(istate))[1:nummod])
		binary_stateTF = unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])

		# Draw from the multivariate normal
		n = length(which(statevec==istate))
		if (n == 0) next
		cat("drawing from multivariate normal...             \r")
		z = matrix(rmvnorm(n, mean=rep(0,nummod), sigma=sigma[,,i]), ncol=nummod)
		# Transform to uniform space
		cat("transforming to uniform space...                \r")
		u = matrix(apply(z, 2, pnorm), ncol=nummod)
		# Transform to count space using marginals
		cat("transform to count space...                     \r")
		isizes = sizes[binary_stateTF]
		iprobs = probs[binary_stateTF]
		for (imod in 1:nummod) {
			mask = statevec==istate
			if (binary_state[imod] == 0) {
				reads[mask,imod] = qzinbinom(u[,imod], w=ws[imod], size=isizes[imod], prob=iprobs[imod])
			} else {
				reads[mask,imod] = qnbinom(u[,imod], size=isizes[imod], prob=iprobs[imod])
			}
		}
		cat("                                                \r")
			
	}

	# Return the output
	out = list(coordinates = coordinates,
				states = statevec,
				reads = reads,
				emissions = emissions,
				transition = transition,
				sigma = sigma,
				use.states = use.states,
				weights = weights
				)
	return(out)

}

simulate.binned2sam <- function(binned.data, chrom.length.file=NULL, file="simulated", fragment.length=1) {
	
	## Assign variables
	names(binned.data) <- binned.data.names
	numbins <- nrow(binned.data)
	numreads <- sum(binned.data$reads)
	flag.labels <- c(0,16)
	chroms <- unique(binned.data$chrom)
	file <- paste(file,"sam", sep=".")
	# Chromosome lengths
	if (!is.null(chrom.length.file)) {
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
	} else {
		chrom.lengths <- rep(NA, length(chroms))
		names(chrom.lengths) <- chroms
		for (ichrom in chroms) {
			iends <- binned.data$end[binned.data$chrom==ichrom]
			chrom.lengths[ichrom] <- iends[length(iends)]
		}
	}

	## Write header
	cat("@PG\tID:chromstar\tVN:0\tCL:simulate.binned2sam()\n", file=file)
	for (ichrom in names(chrom.lengths)) {
		cat(paste("@SQ\tSN:",ichrom,"\tLN:",chrom.lengths[ichrom],"\n", sep=""), file=file, append=T)
	}

	## Make columns for SAM file
	qname <- 1:numreads
	# Go through each chromosome
	icount <- 1
	for (ichrom in chroms) {
		cat(ichrom,"\n")
		imask <- binned.data$chrom==ichrom
		inumreads <- sum(binned.data$reads[imask])
		# Generate SAM columns for chromosome
		cat("\rgenerating reads ...        ")
		samcol <- data.frame(qname = qname[icount:(icount+inumreads-1)],
												flag = sample(flag.labels, inumreads, replace=T),
												rname = rep(ichrom, inumreads),
												pos = rep(0, inumreads),
												mapq = rep(255, inumreads),
												cigar = rep(paste(fragment.length,"M", sep=""), inumreads),
												rnext = rep("*", inumreads),
												pnext = rep(0, inumreads),
												tlen = rep(0, inumreads),
												seq = rep("*", inumreads),
												qual = rep("*", inumreads)
											)
		icount <- icount + inumreads
		samcol$pos <- unlist(apply(binned.data[imask,2:4], 1, function(row) { row<-as.integer(row); as.integer(runif(row[3], row[1]+1, row[2]+1-fragment.length)) } ))
		cat("\rsorting reads ...         ")
		samcol <- samcol[order(samcol$pos),]
		# Removing pos=NA reads, happens when fragment.length is bigger than bin (e.g. in last bin)
		if (any(is.na(samcol$pos))) {
			samcol <- samcol[-which(is.na(samcol$pos)),]
		}
		cat("\rwriting to file ...       ")
		write.table(format(samcol, scientific=F, trim=T), quote=F, sep="\t", row.names=F, col.names=F, append=T, file=file)
		cat("\r                          ")
	}
	cat("\rdone\n")
}

simulate.binned2bed <- function(binned.data, file="simulated", fragment.length=1) {
	
	## Assign variables
	names(binned.data) <- binned.data.names
	numbins <- nrow(binned.data)
	numreads <- sum(binned.data$reads)
	chroms <- unique(binned.data$chrom)
	file <- paste(file,"bed", sep=".")

	## Write header
	# no header

	## Make columns for BED file
	qname <- 1:numreads
	# Go through each chromosome
	icount <- 1
	for (ichrom in chroms) {
		cat(ichrom,"\n")
		imask <- binned.data$chrom==ichrom
		inumreads <- sum(binned.data$reads[imask])
		# Generate BED columns for chromosome
		cat("\rgenerating reads ...        ")
		bedcol <- data.frame(chrom = rep(ichrom, inumreads),
												chromStart = rep(0, inumreads),
												chromEnd = rep(0, inumreads),
												name = rep("N", inumreads),
												score = rep(1000, inumreads),
												strand = sample(c("+","-"), inumreads, replace=T)
											)
		icount <- icount + inumreads
		bedcol$chromStart <- unlist(apply(binned.data[imask,2:4], 1, function(row) { row<-as.integer(row); as.integer(runif(row[3], row[1]+1, row[2]+1-fragment.length)) } ))
		bedcol$chromEnd <- bedcol$chromStart + fragment.length
		cat("\rsorting reads ...         ")
		bedcol <- bedcol[order(bedcol$chromStart),]
		# Removing pos=NA reads, happens when fragment.length is bigger than bin (e.g. in last bin)
		if (any(is.na(bedcol$chromStart))) {
			bedcol <- bedcol[-which(is.na(bedcol$chromStart)),]
		}
		cat("\rwriting to file ...       ")
		write.table(format(bedcol, scientific=F, trim=T), quote=F, sep="\t", row.names=F, col.names=F, append=T, file=file)
		cat("\r                          ")
	}
	cat("\rdone\n")
}
