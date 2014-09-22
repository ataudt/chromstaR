stateorderByTransition <- function(multi.hmm) {

	## Intercept user input
	if (check.multivariate.model(multi.hmm)!=0) {
		cat("Loading multivariate HMM from file ...")
		multi.hmm <- get(load(multi.hmm))
		cat(" done\n")
		if (check.multivariate.model(multi.hmm)!=0) stop("argument 'multi.hmm' expects a multivariate hmm object or a file that contains a multivariate hmm (type ?multi.hmm for help)")
	}

	## Calculate distance matrix
	distances <- matrix(NA, ncol=ncol(multi.hmm$A), nrow=nrow(multi.hmm$A))
	colnames(distances) <- colnames(multi.hmm$A)
	rownames(distances) <- rownames(multi.hmm$A)
	for (irow in 1:nrow(distances)) {
		for (icol in 1:ncol(distances)) {
			distances[irow,icol] <- 2 - (multi.hmm$A[irow,icol] + multi.hmm$A[icol,irow])
# 			distances[irow,icol] <- abs(multi.hmm$A[irow,icol] - multi.hmm$A[icol,irow])
		}
	}

	## Select ordering
	states <- multi.hmm$comb.states
	stateorders <- matrix(NA, ncol=length(states), nrow=length(states))
	total.distance <- rep(0, length(states))
	for (i1 in 1:length(states)) {
		state1 <- states[i1]
		stateorders[i1,1] <- state1
		remaining.states <- setdiff(multi.hmm$comb.states, state1)
		for (i2 in 2:length(states)) {
			next.state <- remaining.states[which.min(distances[as.character(state1),as.character(remaining.states)])]
			stateorders[i1,i2] <- next.state
			next.distance <- min(distances[as.character(state1),as.character(remaining.states)])
			total.distance[i1] <- total.distance[i1] + next.distance
			remaining.states <- setdiff(remaining.states,next.state) 
		}
	}
	stateorder <- stateorders[which.min(total.distance),]

	## Reorder the multi.hmm
	return(stateorder)

}
