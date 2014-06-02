prepare.multivariate = function(modellist, use.states=NULL, num.states=NULL, bernoulli=FALSE) {

	# Function
	initializeTxtProgressBar = function(max) {
		pb = txtProgressBar(min=1, max=max, style=3)
		ipb <<- 1
		setTxtProgressBar(pb, 1)
		return(pb)
	}
	increaseTxtProgressBar = function(pb, ipb) {
		setTxtProgressBar(pb, ipb)
		ipb = ipb+1
		return(ipb)
	}

	nummod = length(modellist)
	numbins = length(modellist[[1]]$reads)

	# Extract the reads
	cat("Extracting reads from modellist...")
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	for (imod in 1:nummod) {
		reads[,imod] = modellist[[imod]]$reads
	}
	maxreads = max(reads)
	cat(" done\n")

	if (bernoulli) {
		# Extract posteriors
		cat("Extracting prob.modified from modellist...")
		prob.modified = matrix(NA, ncol=nummod, nrow=numbins)
		for (imod in 1:nummod) {
			prob.modified[,imod] = modellist[[imod]]$posteriors[,3]
		}
		cat(" done\n")
	} else {
		prob.modified = NULL
	}

	# Extract coordinates and other stuff
	coordinates = modellist[[1]]$coordinates
	distributions = lapply(modellist,"[[","distributions")
	weights = lapply(modellist,"[[","weights")

	# Get the combinatorial states
	cat("Getting combinatorial states...")
	combstates.per.bin = combinatorial.states(modellist)
	comb.states.table = table(combstates.per.bin)
	if (is.null(use.states)) {
		comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	} else {
		comb.states = use.states
	}
	cat(" done\n")

# 	# Calculate a blacklist with bins where one modification has state zeroinflation
# 	bl.states = NULL
# 	for (i1 in 1:nummod) {
# 		bl.states[[i1]] = get.states(modellist[[i1]], aggregate=FALSE)
# 	}
# 	bl.states = matrix(unlist(bl.states), ncol=nummod)
# 	blacklist = rep(FALSE, nrow(bl.states))
# 	for (i1 in 1:nummod) {
# 		blacklist = bl.states[,i1]==-1 | blacklist
# 	}

	# Clean up to reduce memory usage
	remove(modellist)

	## Calculate transition matrix
	cat("Estimating new transition matrix...")
	A.estimated = matrix(0, ncol=2^nummod, nrow=2^nummod)
	colnames(A.estimated) = 1:2^nummod-1
	rownames(A.estimated) = 1:2^nummod-1
	for (i1 in 1:(length(combstates.per.bin)-1)) {
		from = combstates.per.bin[i1] + 1
		to = combstates.per.bin[i1+1] + 1
		A.estimated[from,to] = A.estimated[from,to] + 1
	}
	A.estimated = sweep(A.estimated, 1, rowSums(A.estimated), "/")
	cat(" done\n")
		
	## We pre-compute the z-values for each number of reads
	cat("Computing pre z-matrix...")
	z_per_read = matrix(rep(NA,(maxreads+1)*nummod*2), ncol=nummod*2)
	for (imod in 1:nummod) {
		# Go through unmodified and modified
		for (i1 in 2:3) {
			
			xreads = 0:maxreads
			r = distributions[[imod]][i1,'size']
			p = distributions[[imod]][i1,'prob']
# 			lGamma1plusRplusX = lgamma(1+r+(0:maxreads))
# 			lGammaR = lgamma(r)
# 			lGamma2plusX = lgamma(2+(0:maxreads))
# 			lHyper = log( hyperg_2F1( 1, 1+r+(0:maxreads), 2+(0:maxreads), 1-p ) )
# 			lppowert = (1+(0:maxreads))*log(1-p)
# 			lppowerr = r*log(p)

			if (i1 == 2) {
				# Unmodified with zero inflation
				w = weights[[imod]][1]
				u = pzinbinom(xreads, w, r, p)
# 				u = 1 - exp( log(1-w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )
			} else if (i1 == 3) {
				# Modified
# 				u = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )
				u = pnbinom(xreads, r, p)
			}

			# Check for infinities in u and set them to max value which is not infinity
			qnorm_u = qnorm(u)
			if (NaN %in% qnorm_u) { print("NA detected") }
			testvec = qnorm_u!=Inf
			z_per_read[1:(maxreads+1),(imod - 1) * 2 + (i1-1)] = ifelse(testvec, qnorm_u, max(qnorm_u[testvec]))

		}
	}
	cat(" done\n")

	## Compute the z matrix
	cat("Transfering values into z-matrix...")
# 	pb = initializeTxtProgressBar(nummod*2)
	z_per_bin = matrix(rep(NA,numbins*nummod*2), ncol=nummod*2) # (numbins x 2*nummod) matrix, contains the z-value for each bin and modification enriched and unenriched
	for (imod in 1:nummod) {
		for (i1 in 1:2) {
			z_per_bin[,(imod-1)*2+i1] = z_per_read[reads[,imod]+1,(imod-1)*2+i1]
# 			ipb = increaseTxtProgressBar(pb, ipb)
		}
	}
# 	close(pb)
	cat(" done\n")

	# Clean up to reduce memory usage
	remove(z_per_read)

	## Compute inverse of the correlation matrix
	cat("Computing inverse of correlation matrix...")
# 	pb = initializeTxtProgressBar(length(comb.states))
	z_final = matrix(rep(NA, numbins*nummod), ncol=nummod) # (numbins x nummod) matrix, contains the z-value for each bin and the correct combinatorial state (aka the correct enriched-unenriched combination)
	covarianceMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)))
	correlationMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)))
	correlationMatrixInverse = array(NA, dim=c(nummod,nummod,length(comb.states)))
	determinant = rep(NA, length(comb.states))
	usestateTF = rep(NA,length(comb.states)) # TRUE, FALSE vector for usable states
	for (istate in comb.states) {
		i = which(comb.states==istate)
		
		# Convert istate to binary representation
		binary_state = rev(as.integer(intToBits(istate))[1:nummod])
		binary_stateTF = unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])
		
		# Calculate the correlation matrix
# 		mask = which(combstates.per.bin==istate & !blacklist)
		mask = which(combstates.per.bin==istate)
		z_final[mask,] = z_per_bin[mask, which(binary_stateTF==TRUE)]
		temp = tryCatch({
			covarianceMatrix[,,i] = cov(z_final[mask,])
			correlationMatrix[,,i] = cor(z_final[mask,])
			determinant[i] = det( correlationMatrix[,,i] )
			correlationMatrixInverse[,,i] = solve(correlationMatrix[,,i])
			usestateTF[i] = TRUE
		}, warning = function(war) {
			usestateTF[i] <<- FALSE
			war
		}, error = function(err) {
			usestateTF[i] <<- FALSE
			err
		})
# 		ipb = increaseTxtProgressBar(pb, ipb)
	}
	# Total correlation
	correlationMatrixAllInverse = NA
	temp = tryCatch({
		correlationMatrixAll = cor(z_final)
		determinantAll = det( correlationMatrixAll )
		correlationMatrixAllInverse = solve(correlationMatrixAll)
	}, warning = function(war) {
		war
	}, error = function(err) {
		err
	})
# 	close(pb)
	cat(" done\n")

	# Clean up to reduce memory usage
	remove(z_per_bin)
# 	remove(z_final)

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
		covarianceMatrix2use = covarianceMatrix[,,usestateTF]
		correlationMatrix2use = correlationMatrix[,,usestateTF]
		correlationMatrixInverse2use = correlationMatrixInverse[,,usestateTF]
	} else {
		covarianceMatrix2use = covarianceMatrix[,,usestateTF][,,1:numstates2use]
		correlationMatrix2use = correlationMatrix[,,usestateTF][,,1:numstates2use]
		correlationMatrixInverse2use = correlationMatrixInverse[,,usestateTF][,,1:numstates2use]
	}
	comb.states2use = comb.states[usestateTF][1:numstates2use]
	comb.states.table2use = comb.states.table[as.character(comb.states2use)]
	determinant2use = determinant[usestateTF][1:numstates2use]
	A.estimated2use = A.estimated[as.character(comb.states2use), as.character(comb.states2use)]
	A.estimated2use = sweep(A.estimated2use, 1, rowSums(A.estimated2use), "/") # rescale to rowSums = 1 because of the rows and columns taken out

	# Return parameters
	out = list(coordinates = coordinates,
				reads = reads,
				prob.unmodified = 1-prob.modified,
				numbins = numbins,
				nummod = nummod,
				comb.states = comb.states2use,
				comb.states.table = comb.states.table2use,
				distributions = distributions,
				weights = weights,
				covarianceMatrix = covarianceMatrix2use,
				correlationMatrix = correlationMatrix2use,
				correlationMatrixInverse = correlationMatrixInverse2use,
				determinant = determinant2use,
				correlationMatrixAll = correlationMatrixAll,
				correlationMatrixAllInverse = correlationMatrixAllInverse,
				determinantAll = determinantAll,
				usestateTF = usestateTF,
				numstates2use = numstates2use,
				comb.states.per.bin = combstates.per.bin,
				z = z_final,
				A.estimated = A.estimated2use
	)
	return(out)
}
