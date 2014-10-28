prepare.multivariate = function(modellist, use.states=NULL, num.states=NULL, num.threads=1) {

	nummod = length(modellist)
	numbins = modellist[[1]]$num.bins
	IDs <- unlist(lapply(modellist, "[[", "ID"))

	### Extract the reads
	cat("Extracting reads from modellist...")
	reads = matrix(NA, ncol=nummod, nrow=numbins)
	colnames(reads) <- IDs
	for (imod in 1:nummod) {
		reads[,imod] = modellist[[imod]]$reads
	}
	maxreads = max(reads)
	cat(" done\n")

	### Extract coordinates and other stuff
	coordinates = modellist[[1]]$coordinates
	seqlengths <- modellist[[1]]$seqlengths
	distributions = lapply(modellist,"[[","distributions")
	weights = lapply(modellist,"[[","weights")

	### Get the combinatorial states
	cat("Getting combinatorial states...")
	combstates.per.bin = combinatorial.states(modellist)
	comb.states.table = table(combstates.per.bin)
	if (is.null(use.states)) {
		comb.states = as.numeric(names(sort(comb.states.table, decreasing=TRUE)))
	} else {
		comb.states = use.states
	}
	cat(" done\n")

	# Clean up to reduce memory usage
	remove(modellist)

	## We pre-compute the z-values for each number of reads
	cat("Computing pre z-matrix...")
	z_per_read = matrix(rep(NA,(maxreads+1)*nummod*2), ncol=nummod*2)
	xreads = 0:maxreads
	for (imod in 1:nummod) {
		# Go through unmodified and modified
		for (i1 in 2:3) {
			
			r = distributions[[imod]][i1,'size']
			p = distributions[[imod]][i1,'prob']

			if (i1 == 2) {
				# Unmodified with zero inflation
				w = weights[[imod]][1] / (weights[[imod]][2] + weights[[imod]][1])
				u = pzinbinom(xreads, w, r, p)
			} else if (i1 == 3) {
				# Modified
				u = pnbinom(xreads, r, p)
			}

			# Check for infinities in u and set them to max value which is not infinity
			qnorm_u = qnorm(u)
			if (NaN %in% qnorm_u) { print("NaN detected") }
			testvec = qnorm_u!=Inf
			z_per_read[1:(maxreads+1),(imod - 1) * 2 + (i1-1)] = ifelse(testvec, qnorm_u, max(qnorm_u[testvec]))

		}
	}
	cat(" done\n")

	## Compute the z matrix
	cat("Transfering values into z-matrix...")
	z_per_bin = matrix(rep(NA,numbins*nummod*2), ncol=nummod*2) # (numbins x 2*nummod) matrix, contains the z-value for each bin and modification enriched and unenriched
	for (imod in 1:nummod) {
		for (i1 in 1:2) {
			z_per_bin[,(imod-1)*2+i1] = z_per_read[reads[,imod]+1,(imod-1)*2+i1]
		}
	}

	# Clean up to reduce memory usage
	remove(z_per_read)
	cat(" done\n")

	### Calculate correlation matrix
	cat("Computing inverse of correlation matrix...")
# ptm <- proc.time()
# 	covarianceMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)))
	correlationMatrix = array(NA, dim=c(nummod,nummod,length(comb.states)))
	correlationMatrixInverse = array(NA, dim=c(nummod,nummod,length(comb.states)))
	determinant = rep(NA, length(comb.states))
	usestateTF = rep(NA,length(comb.states)) # TRUE, FALSE vector for usable states

# 	# Setup parallel cluster
# 	library(doParallel)
# 	cl <- makeCluster(num.threads)
# 	registerDoParallel(cl)
# 	result <- foreach (istate = comb.states) %dopar% {
# 		i = which(comb.states==istate)
# 		mask = which(combstates.per.bin==istate)
# 		# Convert istate to binary representation
# 		binary_state = rev(as.integer(intToBits(istate))[1:nummod])
# 		binary_stateTF = unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])
# 		# Subselect z
# 		z_temp <- z_per_bin[mask, which(binary_stateTF==TRUE)]
# 		usestateTF <- NA
# 		determinant <- NA
# 		correlationMatrix = array(NA, dim=c(nummod,nummod))
# 		correlationMatrixInverse = array(NA, dim=c(nummod,nummod))
# 		temp = tryCatch({
# # 			covarianceMatrix = cov(z_temp)
# 			correlationMatrix = cor(z_temp)
# 			determinant = det( correlationMatrix )
# 			correlationMatrixInverse = solve(correlationMatrix)
# 			usestateTF = TRUE
# 		}, warning = function(war) {
# 			usestateTF <<- FALSE
# 			war
# 		}, error = function(err) {
# 			usestateTF <<- FALSE
# 			err
# 		})
# 		return(list(correlationMatrix=correlationMatrix,
# 								determinant=determinant,
# 								correlationMatrixInverse=correlationMatrixInverse,
# 								usestateTF=usestateTF))
# 	}
# 	stopCluster(cl)
# 	# Clean up to reduce memory usage
# 	remove(z_per_bin)
# 
# 	## Transfer from list to array
# 	for (istate in comb.states) {
# 		i = which(comb.states==istate)
# 		correlationMatrix[,,i] = result[[i]]$correlationMatrix
# 		correlationMatrixInverse[,,i] = result[[i]]$correlationMatrixInverse
# 		determinant[i] = result[[i]]$determinant
# 		usestateTF[i] = result[[i]]$usestateTF
# 	}

	## Calculate correlation matrix serial
	for (istate in comb.states) {
		i = which(comb.states==istate)
		mask = which(combstates.per.bin==istate)
		# Convert istate to binary representation
		binary_state = rev(as.integer(intToBits(istate))[1:nummod])
		binary_stateTF = unlist(list(c(TRUE,FALSE), c(FALSE,TRUE))[binary_state+1])
		# Subselect z
		z_temp <- z_per_bin[mask, which(binary_stateTF==TRUE)]
		temp = tryCatch({
# 			covarianceMatrix[,,i] = cov(z_temp)
			correlationMatrix[,,i] = cor(z_temp)
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
	}
	remove(z_per_bin)
	cat(" done\n")
# print(proc.time()-ptm)

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
