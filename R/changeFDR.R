#' Change the false discovery rate of a Hidden Markov Model
#'
#' Adjusts the peak calls of a \code{\link{uniHMM}} or \code{\link{multiHMM}} object with the given false discovery rate (FDR).
#'
#' Peaks are called if the posterior for a state (univariate) or sample (multivariate) are \code{>= 1-FDR}.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} or \code{\link{multiHMM}} object with posteriors.
#' @param FDR False discovery rate.
#' @param separate.zeroinflation Only for \code{\link{uniHMM}} objects: If set to TRUE, state 'zero-inflation' will be treated separately, otherwise it will be merged with state 'unmodified'.
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
#' @examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "euratrans",
#'                       "liver-H3K27me3-BN-male-bio1-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'binned <- binReads(bedfile, assembly=rn4_chrominfo, binsize=1000,
#'                   chromosomes='chr12')
#'## Fit the univariate Hidden Markov Model
#'# !Keep posteriors to change the FDR later!
#'hmm <- callPeaksUnivariate(binned, ID='example_H3K27me3', max.time=60, eps=1,
#'                           keep.posteriors=TRUE)
#'## Compare fits with different FDRs
#'plot(changeFDR(hmm, FDR=0.1), type='histogram') + ylim(0,0.25)
#'plot(hmm, type='histogram') + ylim(0,0.25)
#'plot(changeFDR(hmm, FDR=0.9), type='histogram') + ylim(0,0.25)
#'
changeFDR <- function(model, FDR=0.5, separate.zeroinflation=TRUE, averages=TRUE) {

	FDR <- suppressWarnings( as.numeric(FDR) )
	if (!is.numeric(FDR) | is.na(FDR)) {
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
		ptm <- startTimedMessage("Calculating states from posteriors ...")
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
		stopTimedMessage(ptm)
		## Redo segmentation
		ptm <- startTimedMessage("Making segmentation ...")
		gr <- model$bins
		df <- as.data.frame(gr)
		if (averages==TRUE) {
			red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2average=c('counts','score','posteriors.P.modified.'), columns2drop=c('width','posteriors.P.zero.inflation.','posteriors.P.unmodified.')))
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], mean.counts=red.df[,'mean.counts'], state=red.df[,'state'], score=red.df[,'mean.score'], mean.posterior.modified=red.df[,'mean.posteriors.P.modified.'])
		} else {
			red.df <- suppressMessages(collapseBins(df, column2collapseBy='state', columns2drop=c('width','posteriors.P.zero.inflation.','posteriors.P.unmodified.','posteriors.P.modified.','counts')))
			red.gr <- GRanges(seqnames=red.df[,1], ranges=IRanges(start=red.df[,2], end=red.df[,3]), strand=red.df[,4], state=red.df[,'state'])
		}	
		model$segments <- red.gr
		seqlengths(model$segments) <- seqlengths(model$bins)
		stopTimedMessage(ptm)
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)

	### Multivariate HMM ###
	} else if (is(model,class.multivariate.hmm)) {
		if (is.null(model$bins$posteriors)) stop("Cannot recalculate states because posteriors are missing. Run 'callPeaksMultivariate' again with option 'keep.posteriors' set to TRUE.")
		## Calculate states
		ptm <- startTimedMessage("Calculating states from posteriors ...")
		states <- factor(bin2dec(model$bins$posteriors >= threshold), levels=levels(model$bins$state))
		model$bins$state <- states
		## Combinations
		if (!is.null(model$bins$combination)) {
			mapping <- model$mapping
# 			mapping <- levels(model$bins$combination)
# 			names(mapping) <- levels(model$bins$state)
			model$bins$combination <- factor(mapping[as.character(model$bins$state)], levels=mapping[as.character(levels(model$bins$state))])
		}
		stopTimedMessage(ptm)
		## Redo segmentation
		ptm <- startTimedMessage("Making segmentation ...")
		df <- as.data.frame(model$bins)
		ind.readcols <- grep('^counts', names(df))
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
		stopTimedMessage(ptm)
# 		## Redo weights
# 		model$weights <- table(model$bins$state) / length(model$bins)
	} else {
		stop("Supply either a univariate or multivariate chromstaR model")
	}
	## Return model
	return(model)
	
}

