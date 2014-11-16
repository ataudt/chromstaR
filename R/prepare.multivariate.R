prepare.multivariate = function(modellist, use.states=NULL, num.states=NULL, num.threads=1) {

	nummod = length(modellist)
	numbins = modellist[[1]]$num.bins
	IDs <- unlist(lapply(modellist, "[[", "ID"))

	### Extract the reads
	cat("Extracting reads from modellist...")
	ptm <- proc.time()
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	colnames(reads) <- IDs
	for (imod in 1:nummod) {
		reads[,imod] = modellist[[imod]]$reads
	}
	maxreads = max(reads)
	time <- proc.time() - ptm
	cat(paste0(" ",round(round(time[3], 2),2),"s\n"))

	### Extract coordinates and other stuff
	coordinates = modellist[[1]]$coordinates
	seqlengths <- modellist[[1]]$seqlengths
	distributions = lapply(modellist,"[[","distributions")
	weights = lapply(modellist,"[[","weights")

	### Get the combinatorial states
	cat("Getting combinatorial states...")
	ptm <- proc.time()
	combstates.per.bin = combinatorial.states(modellist)
	comb.states.table = table(combstates.per.bin)
	if (is.null(use.states)) {
		comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	} else {
		comb.states = use.states
	}
	time <- proc.time() - ptm
	cat(paste0(" ",round(time[3], 2),"s\n"))

	# Clean up to reduce memory usage
	remove(modellist)

	## We pre-compute the z-values for each number of reads
	cat("Computing pre z-matrix...")
	ptm <- proc.time()
	z.per.read <- array(NA, dim=c(maxreads+1, nummod, 2))
# 	z_per_read = matrix(rep(NA,(maxreads+1)*nummod*2), ncol=nummod*2)
	xreads = 0:maxreads
	for (imod in 1:nummod) {
		# Go through unmodified and modified
		for (i1 in 2:3) {
			
			size = distributions[[imod]][i1,'size']
			prob = distributions[[imod]][i1,'prob']

			if (i1 == 2) {
				# Unmodified with zero inflation
				w = weights[[imod]][1] / (weights[[imod]][2] + weights[[imod]][1])
				u = pzinbinom(xreads, w, size, prob)
			} else if (i1 == 3) {
				# Modified
				u = pnbinom(xreads, size, prob)
			}

			# Check for infinities in u and set them to max value which is not infinity
			qnorm_u = qnorm(u)
			mask <- qnorm_u==Inf
			qnorm_u[mask] <- qnorm(1-1e-16)
# 			testvec = qnorm_u!=Inf
# 			qnorm_u = ifelse(testvec, qnorm_u, max(qnorm_u[testvec]))
			z.per.read[ , imod, i1-1] <- qnorm_u

		}
	}
	time <- proc.time() - ptm
	cat(paste0(" ",round(time[3], 2),"s\n"))

	## Compute the z matrix
	cat("Transfering values into z-matrix...")
	ptm <- proc.time()
	z.per.bin = array(NA, dim=c(numbins, nummod, 2), dimnames=list(bin=1:numbins, track=IDs, state.labels[2:3]))
	for (imod in 1:nummod) {
		for (i1 in 1:2) {
			z.per.bin[ , imod, i1] <- z.per.read[reads[,imod]+1, imod, i1]
		}
	}

	# Clean up to reduce memory usage
	remove(z.per.read)
	time <- proc.time() - ptm
	cat(paste0(" ",round(time[3], 2),"s\n"))

	### Calculate correlation matrix
	cat("Computing inverse of correlation matrix...")
	ptm <- proc.time()
# 	covarianceMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)))
	correlationMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)), dimnames=list(track=IDs, track=IDs, comb.state=comb.states))
	correlationMatrixInverse = array(NA, dim=c(nummod,nummod,length(comb.states)), dimnames=list(track=IDs, track=IDs, comb.state=comb.states))
	determinant = rep(NA, length(comb.states))
	names(determinant) <- comb.states
	usestateTF = rep(NA,length(comb.states)) # TRUE, FALSE vector for usable states

	## Calculate correlation matrix serial
	for (state in comb.states) {
		istate = which(comb.states==state)
		mask = which(combstates.per.bin==state)
		# Convert state to binary representation
		binary_state = rev(as.integer(intToBits(state))[1:nummod])
		# Subselect z
		z.temp <- matrix(NA, ncol=nummod, nrow=length(mask))
		for (i1 in 1:length(binary_state)) {
			z.temp[,i1] <- z.per.bin[mask, i1, binary_state[i1]+1]
		}
		temp = tryCatch({
# 			covarianceMatrix[,,istate] = cov(z_temp)
			if (nrow(z.temp) > 1) {
				correlationMatrix[,,istate] = cor(z.temp)
				determinant[istate] = det( correlationMatrix[,,istate] )
				correlationMatrixInverse[,,istate] = solve(correlationMatrix[,,istate])
				usestateTF[istate] = TRUE
			} else {
				usestateTF[istate] = FALSE
			}
		}, warning = function(war) {
			usestateTF[istate] <<- FALSE
			war
		}, error = function(err) {
			usestateTF[istate] <<- FALSE
			err
		})
	}
	remove(z.per.bin)
	time <- proc.time() - ptm
	cat(paste0(" ",round(time[3], 2),"s\n"))

	# Determine how many (numstates2use) of the usable (usestateTF) states to use
	ok.numstates = length(which(usestateTF==TRUE))
	max.numstates = length(usestateTF)
	if (is.null(use.states) & is.null(num.states)) {
		numstates2use = ok.numstates
	} else if (!is.null(use.states)) {
		if (ok.numstates < length(use.states)) {
			stop("Cannot use the specified states. The occurrence of the following states is too low: ",paste(use.states[!usestateTF], collapse=" "))
		} else {
			numstates2use = length(use.states)
		}
	} else if (!is.null(num.states)) {
		if (ok.numstates < num.states) {
			if (ok.numstates < max.numstates) {
				warning("Using only ",ok.numstates," states instead of desired ",num.states," states! More is not possible because the occurrence of the following states is too low: ",paste(comb.states[!usestateTF], collapse=" "))
			} else {
				warning("Using only ",ok.numstates," states instead of desired ",num.states," states! More is not possible because more than ",max.numstates," do not exist.")
			}
		}
		numstates2use = min(num.states,ok.numstates)
	}

	# Select only states for which the corMatInv could be calculated (usestateTF==TRUE)
	if (numstates2use == 1) {
# 		covarianceMatrix2use = covarianceMatrix[,,usestateTF]
		correlationMatrix2use = correlationMatrix[,,usestateTF]
		correlationMatrixInverse2use = correlationMatrixInverse[,,usestateTF]
	} else {
# 		covarianceMatrix2use = covarianceMatrix[,,usestateTF][,,1:numstates2use]
		correlationMatrix2use = correlationMatrix[,,usestateTF][,,1:numstates2use]
		correlationMatrixInverse2use = correlationMatrixInverse[,,usestateTF][,,1:numstates2use]
	}
	comb.states2use = comb.states[usestateTF][1:numstates2use]
	comb.states.table2use = comb.states.table[as.character(comb.states2use)]
	determinant2use = determinant[usestateTF][1:numstates2use]




# 	## Calculate multivariate densities for each state
# 	cat("Calculating multivariate densities...")
# 	densities <- matrix(1, ncol=numstates2use, nrow=numbins, dimnames=list(bin=1:numbins, comb.state=comb.states2use))
# 	for (state in comb.states2use) {
# 		i <- which(state==comb.states2use)
# 		binary_state = rev(as.integer(intToBits(state))[1:nummod])
# 		product <- 1
# 		for (imod in 1:nummod) {
# 			cat("track ",imod,"\n")
# 			if (binary_state[imod]==0) {
# 				exponent <- -0.5 * apply( ( z.per.bin[ , , 'unmodified'] %*% (correlationMatrixInverse2use[ , , i] - diag(nummod)) ) * z.per.bin[ , , 'unmodified'], 1, sum)
# 				size <- distributions[[imod]]['unmodified','size']
# 				prob <- distributions[[imod]]['unmodified','prob']
# 				w = weights[[imod]][1] / (weights[[imod]][2] + weights[[imod]][1])
# 				product <- product * dzinbinom(reads[,imod], w, size, prob)
# 				densities[,i] <- product * determinant2use[i]^(-0.5) * exp( exponent )
# 			} else if (binary_state[imod]==1) {
# 				exponent <- -0.5 * apply( ( z.per.bin[ , , 'modified'] %*% (correlationMatrixInverse2use[ , , i] - diag(nummod)) ) * z.per.bin[ , , 'modified'], 1, sum)
# 				size <- distributions[[imod]]['modified','size']
# 				prob <- distributions[[imod]]['modified','prob']
# 				product <- product * dnbinom(reads[,imod], size, prob)
# 				densities[,i] <- product * determinant2use[i]^(-0.5) * exp( exponent )
# 			}
# 		}
# 	}
# 	# Check if densities are > 1
# 	if (any(densities>1)) stop("Densities > 1")
# 	if (any(densities<0)) stop("Densities < 0")
# 	densities[densities>1] <- 1
# 	densities[densities<0] <- 0
# 	# Check if densities are 0 everywhere in some bins
# 	check <- which(apply(densities, 1, sum) == 0)
# 	if (length(check)>0) {
# 		if (check[1]==1) {
# 			densities[1,] <- rep(1e-10, ncol(densities))
# 			check <- check[-1]
# 		}
# 		for (icheck in check) {
# 			densities[icheck,] <- densities[icheck-1,]
# 		}
# 	}
# 	cat(" done\n")



	# Return parameters
	out = list(IDs = IDs,
				coordinates = coordinates,
				seqlengths = seqlengths,
				reads = reads,
				numbins = numbins,
				nummod = nummod,
				comb.states = comb.states2use,
				comb.states.table = comb.states.table2use,
				distributions = distributions,
				weights = weights,
# 				covarianceMatrix = covarianceMatrix2use,
				correlationMatrix = correlationMatrix2use,
				correlationMatrixInverse = correlationMatrixInverse2use,
				determinant = determinant2use,
				usestateTF = usestateTF,
				numstates2use = numstates2use,
				comb.states.per.bin = combstates.per.bin
	)
	return(out)
}
