get.states <- function(model, threshold=0.5, separate.zeroinflation=FALSE, control=FALSE) {
	states <- rep(NA,model$num.bins)
	if (is(model,class.chromstar.univariate)) {
		if (control) {
			states <- ifelse(model$posteriors[,2]>threshold, 2, 1)
			states <- factor(state.labels[1:2], levels=states.labels[1:2])[states]
		} else if (separate.zeroinflation) {
			states[ model$posteriors[,3]<=threshold & model$posteriors[,2]<=model$posteriors[,1] ] <- 1
			states[ model$posteriors[,3]<=threshold & model$posteriors[,2]>=model$posteriors[,1] ] <- 2
			states[ model$posteriors[,3]>threshold ] <- 3
			states <- factor(state.labels, levels=state.labels)[states]
		} else {
			states <- ifelse(model$posteriors[,3]>threshold, 2, 1)
			states <- factor(state.labels[2:3], levels=states.labels[2:3])[states]
		}
		return(states)
	} else if (is(model,class.chromstar.multivariate)) {
		stop("Multivariate get.states() not implemented")
	} else {
		stop("Supply either a univariate or multivariate chromstar model")
	}
	
}

