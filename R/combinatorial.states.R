combinatorial.states = function(modellist, binary=FALSE) {
	
# 	# Combine the input
# 	models = c(list(...), modellist) # This will result in excessive memory usage if the list elements are large, because they will be copied to a new location
	nummod = length(modellist)
	numbins = length(modellist[[1]]$bins)

	# Get the univariate states (zero inflation = 0, unmodified = 0, modified = 1) from the modellist
	binary_statesmatrix = matrix(rep(NA,numbins*nummod), ncol=nummod)
	for (imod in 1:nummod) {
		binary_statesmatrix[,imod] = c(FALSE,FALSE,TRUE)[modellist[[imod]]$bins$state] # F,F,T corresponds to levels 'zero-inflation','unmodified','modified'
	}

	if (binary == TRUE) {
		return(binary_statesmatrix)
	} else {
	
		# Transform binary to decimal
		decimal_states = rep(0,numbins)
		for (imod in 1:nummod) {
			decimal_states = decimal_states + 2^(nummod-imod) * binary_statesmatrix[,imod]
		}

		return(decimal_states)

	}

}
