#' Change the false discovery rate of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} object with the given false discovery rate (FDR).
#'
#' Peaks are called if the posterior for a state (univariate) or sample (multivariate) are \code{>= 1-FDR}.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} object with posteriors.
#' @param FDR False discovery rate.
#' @param separate.zeroinflation Only for \code{\link{chromstaR_univariateHMM}} objects: If set to TRUE, state 'zero-inflation' will be treated separately, otherwise it will be merged with state 'unmodified'.
#' @return The input object is returned with adjusted peak calls.
#' @examples
#'## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue
#'bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz",
#'                       package="chromstaR")
#'## Bin the BED file into bin size 1000bp
#'binned.data <- bed2binned(bedfile, assembly='hg19', binsize=1000,
#'                          save.as.RData=FALSE)
#'## Fit the univariate Hidden Markov Model (and keep posteriors!)
#'hmm <- callPeaksUnivariate(binned.data, ID='example_H3K36me3', max.time=60,
#'                           FDR=0.5, keep.posteriors=TRUE)
#'## Adjust the FDR
#'hmm.adjusted <- changeFDR(hmm, FDR=0.01)
#'## Compare the state calls
#'table(hmm$bins$state)
#'table(hmm.adjusted$bins$state)
#' @export
changeFDR <- function(model, FDR=0.5, separate.zeroinflation=TRUE) {

	## Check if posteriors are present
	if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksUnivariate' again with option 'keep.posteriors' set to TRUE.")

	## Get the states
	threshold <- 1-FDR
	model$FDR <- FDR

	### Univariate HMM ###
	if (is(model,class.univariate.hmm)) {
		## Calculate states
		message("Calculating states from posteriors ...", appendLF=F)
		ptm <- proc.time()
		states <- rep(NA,length(model$bins))
		if (separate.zeroinflation) {
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
		message(" ",round(time[3],2),"s")
		## Redo segmentation
		message("Making segmentation ...", appendLF=F)
		ptm <- proc.time()
		gr <- model$bins
		red.gr.list <- GRangesList()
		for (state in state.labels) {
			red.gr <- GenomicRanges::reduce(gr[states==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(state.labels)),length(red.gr))
			if (length(red.gr)>0) {
				red.gr.list[[length(red.gr.list)+1]] <- red.gr
			}
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)

	### Multivariate HMM ###
	} else if (is(model,class.multivariate.hmm)) {
		message("Calculating states from posteriors ...", appendLF=F)
		ptm <- proc.time()
		if (is.null(FDR)) {
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
		message(" ",round(time[3],2),"s")
		## Redo segmentation
		message("Making segmentation ...", appendLF=F)
		ptm <- proc.time()
		gr <- model$bins
		red.gr.list <- GRangesList()
		for (state in levels(model$bins$state)) {
			red.gr <- GenomicRanges::reduce(gr[gr$state==state])
			mcols(red.gr)$state <- rep(factor(state, levels=levels(gr$state)),length(red.gr))
			if (length(red.gr)>0) {
				red.gr.list[[length(red.gr.list)+1]] <- red.gr
			}
		}
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm
		message(" ",round(time[3],2),"s")
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)
	} else {
		stop("Supply either a univariate or multivariate chromstaR model")
	}
	## Return model
	return(model)
	
}

