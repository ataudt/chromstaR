change.FDR <- function(model, false.discovery.rate=0.5, separate.zeroinflation=TRUE, control=FALSE) {

	## Check if posteriors are present
	if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'call.peaks.univariate' again with option 'keep.posteriors' set to TRUE.")

	## Get the states
	threshold <- 1-false.discovery.rate
	model$FDR <- false.discovery.rate
	if (is(model,class.univariate.hmm)) {
		## Calculate states
		cat("Calculating states from posteriors ...")
		ptm <- proc.time()
		states <- rep(NA,length(model$bins))
		if (control) {
			states <- ifelse(model$bins$posteriors[,2]>threshold, 2, 1)
			states <- factor(state.labels[1:2], levels=states.labels[1:2])[states]
		} else if (separate.zeroinflation) {
			states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]<=model$bins$posteriors[,1] ] <- 1
			states[ model$bins$posteriors[,3]<threshold & model$bins$posteriors[,2]>model$bins$posteriors[,1] ] <- 2
			states[ model$bins$posteriors[,3]>=threshold ] <- 3
			states <- state.labels[states]
		} else {
			states <- ifelse(model$bins$posteriors[,3]>=threshold, 2, 1)
			states <- state.labels[2:3][states]
		}
		model$bins$state <- states
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		## Redo segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- model$bins
		red.gr.list <- GRangesList()
		for (state in state.labels) {
			red.gr <- GenomicRanges::reduce(gr[states==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)
	} else if (is(model,class.multivariate.hmm)) {
		cat("Calculating states from posteriors ...")
		ptm <- proc.time()
		if (is.null(false.discovery.rate)) {
			states <- factor(levels(model$bins$state)[apply(model$bins$posteriors, 1, which.max)], levels=levels(model$bins$state))
		} else {
			post.per.track <- matrix(0, ncol=ncol(model$bins$reads), nrow=nrow(model$bins$reads))
			binstates <- dec2bin(levels(model$bins$state), ndigits=ncol(model$bins$reads))
			for (icol in 1:ncol(post.per.track)) {
				binstate.matrix <- matrix(rep(binstates[icol,], nrow(model$bins$reads)), nrow=nrow(model$bins$reads), byrow=T)
				post.per.track <- post.per.track + binstate.matrix * model$bins$posteriors[,icol]
			}
			states <- factor(bin2dec(post.per.track >= threshold), levels=levels(model$bins$state))
# 			max.index <- apply(model$bins$posteriors, 1, which.max)
# 			max.states <- factor(levels(model$bins$state)[max.index], levels=levels(model$bins$state))
# 			above.threshold <- apply(model$bins$posteriors, 1, function(x) { any(x>=threshold) })
# 			states <- factor(rep(0, length(model$bins)), levels=unique(c(0,levels(model$bins$state))))
# 			states[above.threshold] <- max.states[above.threshold]
		}
		model$bins$state <- states
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
		## Redo segmentation
		cat("Making segmentation ...")
		ptm <- proc.time()
		gr <- model$bins
		red.gr.list <- GRangesList()
		for (state in levels(model$bins$state)) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm
		cat(paste0(" ",round(time[3],2),"s\n"))
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)
	} else {
		stop("Supply either a univariate or multivariate chromstar model")
	}
	## Return model
	return(model)
	
}

