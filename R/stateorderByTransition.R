stateorderByTransition <- function(transitionProbs) {

    ## Calculate distance matrix
    distances <- matrix(NA, ncol=ncol(transitionProbs), nrow=nrow(transitionProbs))
    colnames(distances) <- colnames(transitionProbs)
    rownames(distances) <- rownames(transitionProbs)
    for (irow in 1:nrow(distances)) {
        for (icol in 1:ncol(distances)) {
            distances[irow,icol] <- 2 - (transitionProbs[irow,icol] + transitionProbs[icol,irow])
#             distances[irow,icol] <- abs(transitionProbs[irow,icol] - transitionProbs[icol,irow])
        }
    }

    ## Select ordering
    comb.states <- colnames(transitionProbs)
    stateorders <- matrix(NA, ncol=length(comb.states), nrow=length(comb.states))
    total.distance <- rep(0, length(comb.states))
    for (i1 in 1:length(comb.states)) {
        state1 <- comb.states[i1]
        stateorders[i1,1] <- state1
        remaining.comb.states <- setdiff(comb.states, state1)
        for (i2 in 2:length(comb.states)) {
            next.state <- remaining.comb.states[which.min(distances[as.character(state1),as.character(remaining.comb.states)])]
            stateorders[i1,i2] <- next.state
            next.distance <- min(distances[as.character(state1),as.character(remaining.comb.states)])
            total.distance[i1] <- total.distance[i1] + next.distance
            remaining.comb.states <- setdiff(remaining.comb.states,next.state) 
        }
    }
    stateorder <- stateorders[which.min(total.distance),]
    return(stateorder)

}
