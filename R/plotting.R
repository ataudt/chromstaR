#' @import ggplot2
#' @import reshape2
NULL

# =================================================================
# Define plotting methods for the generic
# =================================================================
#' Plotting function for saved \pkg{\link{chromstaR}} objects
#'
#' Convenience function that loads and plots a \pkg{\link{chromstaR}} object in one step.
#'
#' @param x A filename that contains either \code{\link{binned.data}}, a \code{\link{chromstaR_univariateHMM}} or a \code{\link{chromstaR_multivariateHMM}}.
#' @param ... Additional arguments.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot character
#' @export
plot.character <- function(x, ...) {
	x <- get(load(x))
	plot(x, ...)
}

#' Plotting function for binned read counts
#'
#' Make plots for binned read counts from \code{\link{binned.data}}
#'
#' @param x A \code{\link{GRanges}} object with binned read counts.
#' @inheritParams plotBinnedDataHistogram
#' @param ... Additional arguments not implemented.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot GRanges
#' @export
plot.GRanges <- function(x, chromosomes=NULL, start=NULL, end=NULL, ...) {
	plotBinnedDataHistogram(x, chromosomes=NULL, start=NULL, end=NULL, ...)
}

#' Plotting function for \code{\link{chromstaR_univariateHMM}} objects
#'
#' Make different types of plots for \code{\link{chromstaR_univariateHMM}} objects.
#'
#' @param x A \code{\link{chromstaR_univariateHMM}} object.
#' @param type Type of the plot, one of \code{c('histogram', 'karyogram', 'boxplot', 'normalTransformation')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{A karyogram with binned read counts and peak calls. This uses the \pkg{\link{ggbio}} package and is very slow!}
#'   \item{\code{boxplot}}{A boxplot of read counts for the different states.}
#'   \item{\code{normalTransformation}}{A histogram of transformed read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' \describe{
#'   \item{\code{state}}{Plot the \code{histogram}, \code{boxplot} or \code{normalTransformation} only for the specified state. One of \code{c('unmodified','modified')}.}
#'   \item{\code{chromosomes,start,end}}{Plot the \code{histogram} only for the specified chromosomes, start and end position.}
#' }
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot chromstaR_univariateHMM
#' @export
plot.chromstaR_univariateHMM <- function(x, type='histogram', ...) {
	
	if (type == 'histogram' | type==1) {
		plotUnivariateHistogram(x, ...)
	} else if (type == 'karyogram' | type==2) {
		plotUnivariateKaryogram(x)
	} else if (type == 'boxplot' | type==3) {
		plotUnivariateBoxplot(x, ...)
	} else if (type == 'normalTransformation' | type==4) {
		plotUnivariateNormalTransformation(x, ...)
	}

}

#' Plotting function for \code{\link{chromstaR_multivariateHMM}} objects
#'
#' Make different types of plots for \code{\link{chromstaR_multivariateHMM}} objects.
#'
#' @param x A \code{\link{chromstaR_multivariateHMM}} object.
#' @param type Type of the plot, one of \code{c('transitionMatrix','histograms')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{transitionMatrix}}{A heatmap with entries of the transition matrix.}
#'   \item{\code{histograms}}{Fitted histograms of all underlying univariate distributions.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot chromstaR_multivariateHMM
#' @export
plot.chromstaR_multivariateHMM <- function(x, type='transitionMatrix', ...) {

	if (type == 'transitionMatrix' | type==1) {
		plotMultivariateTransition(x, ...)
	} else if (type == 'histograms' | type==2) {
		plotMultivariateHistograms(x, ...)
	}

}

# ============================================================
# Plot a read histogram
# ============================================================
#' Plot a histogram of binned read counts
#'
#' Plot a histogram of binned read counts from \code{\link{binned.data}}
#'
#' @param binned.data A \code{\link{binned.data}} object containing binned read counts in meta-column 'reads'.
#' @param chromosomes,start,end Plot the histogram only for the specified chromosomes, start and end position.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotBinnedDataHistogram <- function(binned.data, chromosomes=NULL, start=NULL, end=NULL) {

	# -----------------------------------------
	# Get right x limit
	get_rightxlim <- function(histdata, reads) {
		if (!any(!reads==0)) {
			return(rightxlim=10)
		}
		rightxlim1 <- median(reads[reads>0])*7
		breaks <- histdata$breaks[1:length(histdata$counts)]
		counts <- histdata$counts
		rightxlim2 <- breaks[counts<=5 & breaks>median(reads)*2][1]
		rightxlim <- min(c(rightxlim1,rightxlim2), na.rm=TRUE)
		return(rightxlim)
	}

	# Select the rows to plot
	selectmask <- rep(TRUE,length(binned.data))
	numchrom <- length(table(seqnames(binned.data)))
	if (!is.null(chromosomes)) {
		if (any(! chromosomes %in% levels(seqnames(binned.data)))) {
			stop(chromosomes[! chromosomes %in% levels(seqnames(binned.data))]," can't be found in the binned data.")
		}
		selectchrom <- seqnames(binned.data) %in% chromosomes
		selectmask <- selectmask & selectchrom
		numchrom <- 1
	}
	if (numchrom == 1) {
		if (!is.null(start)) {
			selectstart <- start(ranges(binned.data)) >= start
			selectmask <- selectmask & selectstart
		}
		if (!is.null(end)) {
			selectend <- end(ranges(binned.data)) <= end
			selectmask <- selectmask & selectend
		}
	}
	if (length(which(selectmask)) != length(binned.data$reads)) {
		reads <- binned.data$reads[selectmask]
	} else {
		reads <- binned.data$reads
	}

	# Find the x limits
	breaks <- max(reads)
	if (max(reads)==0) { breaks <- 1 }
	histdata <- hist(reads, right=FALSE, breaks=breaks, plot=FALSE)
	rightxlim <- get_rightxlim(histdata, reads)

	# Plot the histogram
	ggplt <- ggplot(data.frame(reads)) + geom_histogram(aes_string(x='reads', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")
	return(ggplt)

}

# =================================================================
# Plot a read histogram with univariate fits for a multivariate HMM
# =================================================================
plotMultivariateHistograms <- function(multi.hmm) {

	## Make fake uni.hmm and plot
	ggplts <- list()
	for (i1 in 1:length(multi.hmm$IDs)) {
		uni.hmm <- list()
		uni.hmm$ID <- multi.hmm$IDs[i1]
		uni.hmm$bins <- multi.hmm$bins
		uni.hmm$bins$state <- NULL
		uni.hmm$bins$reads <- multi.hmm$bins$reads[,i1]
		uni.hmm$weights <- multi.hmm$weights.univariate[[i1]]
		uni.hmm$distributions <- multi.hmm$distributions[[i1]]
		class(uni.hmm) <- class.univariate.hmm
		ggplts[[i1]] <- plotUnivariateHistogram(uni.hmm)
	}
	
	return(ggplts)

}

# ============================================================
# Plot a read histogram with univariate fits
# ============================================================
#' Plot a histogram of binned read counts with fitted mixture distribution
#'
#' Plot a histogram of binned read counts from with fitted mixture distributions from a \code{\link{chromstaR_univariateHMM}} object.
#'
#' @param model A \code{\link{chromstaR_univariateHMM}} object.
#' @param state Plot the histogram only for the specified state. One of \code{c('unmodified','modified')}.
#' @param chromosomes,start,end Plot the histogram only for the specified chromosomes, start and end position.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotUnivariateHistogram <- function(model, state=NULL, chromosomes=NULL, start=NULL, end=NULL) {

	## Check user input
	if (check.univariate.model(model)!=0) {
		model <- get(load(model))
		if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate HMM or a file that contains a univariate HMM")
	}

	# -----------------------------------------
	# Get right x limit
	get_rightxlim <- function(histdata, reads) {
		if (!any(!reads==0)) {
			return(rightxlim=10)
		}
		rightxlim1 <- median(reads[reads>0])*7
		breaks <- histdata$breaks[1:length(histdata$counts)]
		counts <- histdata$counts
		rightxlim2 <- breaks[counts<=5 & breaks>median(reads)][1]
		rightxlim <- min(c(rightxlim1,rightxlim2), na.rm=TRUE)
		return(rightxlim)
	}

	# Select the rows to plot
	selectmask <- rep(TRUE,length(model$bins))
	if (!is.null(chromosomes)) {
		if (any(! chromosomes %in% levels(seqnames(model$bins)))) {
			stop(chromosomes[! chromosomes %in% levels(seqnames(model$bins))]," can't be found in the binned data.")
		}
		selectchrom <- as.logical(seqnames(model$bins) %in% chromosomes)
		selectmask <- selectmask & selectchrom
		if (!is.null(start)) {
			selectstart <- as.logical(start(ranges(model$bins)) >= start)
			selectmask <- selectmask & selectstart
		}
		if (!is.null(end)) {
			selectend <- as.logical(end(ranges(model$bins)) <= end)
			selectmask <- selectmask & selectend
		}
	}
	if (!is.null(state)) {
		selectmask <- selectmask & model$bins$state==state
	}
	if (length(which(selectmask)) != length(model$bins$reads)) {
		reads <- model$bins$reads[selectmask]
		states <- model$bins$state[selectmask]
		weights <- rep(NA, 3)
		weights[1] <- length(which(states=="zero-inflation"))
		weights[2] <- length(which(states=="unmodified"))
		weights[3] <- length(which(states=="modified"))
		weights <- weights / length(states)
	} else {
		reads <- model$bins$reads
		weights <- model$weights
	}

	# Find the x limits
	breaks <- max(reads)
	if (max(reads)==0) { breaks <- 1 }
	histdata <- hist(reads, right=FALSE, breaks=breaks, plot=FALSE)
	rightxlim <- get_rightxlim(histdata, reads)

	# Plot the histogram
	ggplt <- ggplot(data.frame(reads)) + geom_histogram(aes_string(x='reads', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")

	### Add fits to the histogram
	x <- 0:max(reads)
	distributions <- data.frame(x)

	# Unmodified
	w <- weights[1]/(weights[2]+weights[1]) # weight for the zero-inflation
	if (is.nan(w)) { w <- 0 }
	distributions$unmodified <- (1-weights[3]) * dzinbinom(x, w, model$distributions[2,'size'], model$distributions[2,'prob'])
	# Modified
	distributions$modified <- weights[3] * dnbinom(x, model$distributions[3,'size'], model$distributions[3,'prob'])
	# Total
	distributions$total <- distributions$unmodified + distributions$modified
	# Convert to long format
	df <- melt(distributions, id.vars='x', variable.name='state', value.name='y')

	# Make legend
	lmeans <- round(model$distributions[,'mu'], 2)[-1]
	lvars <- round(model$distributions[,'variance'], 2)[-1]
	lweights <- round(c(1-weights[3], weights[3]), 2)
	legend <- paste0(c('unmodified','modified'), ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
	legend <- c(legend, paste0('total, mean(data)=', round(mean(reads),2), ', var(data)=', round(var(reads),2)))
	ggplt <- ggplt + ggtitle(model$ID)

	### Plot the distributions
	if (is.null(state)) {
		ggplt <- ggplt + geom_line(data=df, aes_string(x='x', y='y', col='state'))
		ggplt <- ggplt + scale_color_manual(name="components", values=state.colors[c('unmodified','modified','total')], labels=legend) + theme(legend.justification=c(1,1), legend.position=c(1,1))
	} else {
		if (state=="unmodified") {
			ggplt <- ggplt + geom_line(data=df[df$state=='unmodified',], aes_string(x='x', y='y', col='state'))
			ggplt <- ggplt + scale_color_manual(name="components", values=state.colors[c('unmodified')], labels=legend[1]) + theme(legend.justification=c(1,1), legend.position=c(1,1))
		}
		if (state=="modified") {
			ggplt <- ggplt + geom_line(data=df[df$state=='modified',], aes_string(x='x', y='y', col='state'))
			ggplt <- ggplt + scale_color_manual(name="components", values=state.colors[c('modified')], labels=legend[2]) + theme(legend.justification=c(1,1), legend.position=c(1,1))
		}
	}
		
	return(ggplt)

}

# ============================================================
# Plot a karyogram with reads and univariate calls
# ============================================================
#' Plot a karyogram with read counts and univariate peak calls
#'
#' Plot a karyogram with read counts and peak calls from a \code{\link{chromstaR_univariateHMM}} object.
#'
#' @author Aaron Taudt
#' @param mode A \code{\link{chromstaR_univariateHMM}} object or file that contains such an object.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @export
plotUnivariateKaryogram <- function(model) {

	## Check user input
	if (check.univariate.model(model)!=0) {
		model <- get(load(model))
		if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate HMM or a file that contains a univariate HMM")
	}

	# Plot the peaks
	ggplt <- ggbio::autoplot(model$segments[model$segments$state=='modified'], layout='karyogram', color=state.colors['modified'])

	# Plot the read counts
	ggplt <- ggplt + ggbio::layout_karyogram(model$bins, aes_string(x='start', y='reads'), ylim=c(10,40), geom='line', color=state.colors['reads'])
	ggplt <- ggplt + theme(axis.text.y=element_blank(), panel.background=element_blank())

	return(ggplt)
}

# ============================================================
# Plot a read histogram in normal space of the given state
# ============================================================
plotUnivariateNormalTransformation <- function(model, state='unmodified') {

	## Plot settings
	cols <- state.colors[c("unmodified","modified")]

	## Transform the reads
	df <- as.data.frame(model$bins)[,c('reads','state')]
	# Transform to uniform space
	mask <- df$state=='modified'
	df$ureads[!mask] <- pzinbinom(df$reads[!mask], model$weights[1], model$distributions[2,'size'], model$distributions[2,'prob'])
	df$ureads[mask] <- pnbinom(df$reads[mask], model$distributions[3,'size'], model$distributions[3,'prob'])
	# Transform to normal space
	df$nreads <- qnorm(df$ureads)

	## Make the plots
	subset <- df$nreads[df$state==state]
	breaks <- c(-Inf,sort(as.numeric(names(table(subset)))))
	x <- seq(-4,4,0.1)
	title <- paste0("Transformed emission density for state ",state)
	ggplt <- ggplot() + geom_histogram(data=data.frame(ureads=subset), aes_string(x='ureads', y='..density..'), breaks=breaks, right=TRUE, col='black', fill=cols[state]) + theme_bw() + geom_line(data=data.frame(x=x, y=dnorm(x, mean=0, sd=1)), aes_string(x='x', y='y')) + xlab("transformed reads") + labs(title=title)
	return(ggplt)

}

# ============================================================
# Plot a boxplot of the univariate calls
# ============================================================
plotUnivariateBoxplot <- function(model) {

	## Boxplot
	df <- as.data.frame(model$bins)[,c('state','reads')]
	ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes_string(x='state', y='reads', fill='state')) + scale_fill_manual(values=state.colors)
	return(ggplt)

}

# ============================================================
# Plot a heat map of the transition probabilities
# ============================================================
plotMultivariateTransition <- function(multi.hmm) {

	A <- melt(multi.hmm$transitionProbs, varnames=c('from','to'), value.name='prob')
	A$from <- factor(A$from, levels=stateorderByTransition(multi.hmm))
	A$to <- factor(A$to, levels=stateorderByTransition(multi.hmm))
	ggplt <- ggplot(data=A) + geom_tile(aes_string(x='to', y='from', fill='prob')) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low="white", high="blue")

	return(ggplt)

}

# ============================================================
# Plot a karyogram with combinatorial entropy
# ============================================================
#' Plot combinatorial entropy
#'
#' Plot a karyogram with combinatorial entropy generated by \code{\link{combinatorialEntropy}}.
#'
#' @author Aaron Taudt
#' @param gr A \code{\link{GRanges}} object generated by \code{\link{combinatorialEntropy}}
#' @param genome The genome assembly to use. Type \code{rtracklayer::ucscGenomes()$db} for a list of available genomes.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @export
plotCombinatorialEntropy <- function(gr, genome) {

	## Get cytoband info
	ideogram <- biovizBase::getIdeogram(genome=genome, cytoband=TRUE)
	cytoband <- keepSeqlevels(ideogram, mixedsort(seqlevels(ideogram)))
	
	# Reshape GRanges to get nice plotting
	gr.low <- resize(gr, width=1, fix='start')
	gr.high <- resize(gr, width=1, fix='end')
	gr <- sort(c(gr.low, gr.high))

	## Make the karyogram
	ggcytoband <- ggbio::autoplot(cytoband, layout='karyogram', cytoband=TRUE)
	ggplt <- ggcytoband + ggbio::layout_karyogram(gr, aes_string(x='start', y='entropy'), ylim=c(10,40), geom='line', color='red')
	ggplt <- ggplt + theme(axis.text.y=element_blank(), panel.background=element_blank())

	return(ggplt)
}
