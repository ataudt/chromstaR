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
#' @param averages Whether or not averaged posteriors should appear in the segmentation.
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
changeFDR <- function(model, FDR=0.5, separate.zeroinflation=TRUE, averages=TRUE) {

	if (!is.numeric(FDR)) {
		warning("Not changing FDR because given FDR is not numeric.")
		return(model)
	} else if (FDR < 0 | FDR > 1) {
		warning("FDR has to be inside the interval [0,1]. Nothing done.")
		return(model)
	}

	## FDR threshold
	threshold <- 1-FDR
	model$FDR <- FDR

	### Univariate HMM ###
	if (is(model,class.univariate.hmm)) {
		if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksUnivariate' again with option 'keep.posteriors' set to TRUE.")
		## Calculate states
		message("Calculating states from posteriors ...", appendLF=F); ptm <- proc.time()
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
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		## Redo segmentation
		message("Making segmentation ...", appendLF=F); ptm <- proc.time()
		gr <- model$bins
		df <- as.data.frame(gr)
		if (averages==TRUE) {
			red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('reads','score','posteriors.P.modified.'), columns2drop=c('width','posteriors.P.zero.inflation.','posteriors.P.unmodified.')))
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], mean.reads=red.df[,'mean.reads'], state=red.df[,'state'], score=red.df[,'mean.score'], mean.posterior.modified=red.df[,'mean.posteriors.P.modified.'])
		} else {
			red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width','posteriors.P.zero.inflation.','posteriors.P.unmodified.','posteriors.P.modified.','reads')))
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'])
		}	
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)

	### Multivariate HMM ###
	} else if (is(model,class.multivariate.hmm)) {
		if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksMultivariate' again with option 'keep.posteriors' set to TRUE.")
		## Calculate states
		message("Calculating states from posteriors ...", appendLF=F); ptm <- proc.time()
		states <- factor(bin2dec(model$bins$posteriors >= threshold), levels=levels(model$bins$state))
		model$bins$state <- states
		## Combinations
		if (!is.null(model$bins$combination)) {
			mapping <- model$mapping
# 			mapping <- levels(model$bins$combination)
# 			names(mapping) <- levels(model$bins$state)
			model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		## Redo segmentation
		message("Making segmentation ...", appendLF=F); ptm <- proc.time()
		df <- as.data.frame(model$bins)
		ind.readcols <- grep('^reads', names(df))
		ind.postcols <- grep('^posteriors', names(df))
		ind.widthcol <- grep('width', names(df))
		ind.scorecol <- grep('score', names(df))
		if (averages==TRUE) {
			red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c(ind.postcols,ind.scorecol), columns2drop=c(ind.readcols, ind.widthcol)))
			mean.posteriors <- matrix(unlist(red.df[,grepl('^mean.posteriors',names(red.df))]), ncol=length(model$IDs))
			colnames(mean.posteriors) <- model$IDs
			mean.score <- red.df[,grepl('^mean.score', names(red.df))]
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'], score=mean.score)
			red.gr$mean.posteriors <- mean.posteriors
			if (!is.null(model$bins$combination)) {
				red.gr$combination <- red.df[,'combination']
			}
		} else {
			red.df <- suppressMessages(collapseBins(df[,-c(4, ind.readcols, ind.postcols)], column2collapseBy='state'))
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'])
			if (!is.null(model$bins$combination)) {
				red.gr$combination <- red.df[,'combination']
			}
		}
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)
	} else {
		stop("Supply either a univariate or multivariate chromstaR model")
	}
	## Return model
	return(model)
	
}

