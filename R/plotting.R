# =================================================================
# Define plotting methods for the generic
# =================================================================
plot.chromstaR_univariateHMM <- function(x, type='histogram', ...) {
	
	if (type == 'histogram' | type==1) {
		plotUnivariateHistogram(x, ...)
	} else if (type == 'boxplot' | type==2) {
		plotUnivariateBoxplot(x, ...)
	} else if (type == 'normalTransformation' | type==3) {
		plotUnivariateNormalTransformation(x, ...)
	}

}

plot.chromstaR_multivariateHMM <- function(x, type='transitionMatrix', ...) {

	if (type == 'transitionMatrix' | type==1) {
		plotMultivariateTransition(x, ...)
	} else if (type == 'histograms' | type==2) {
		plotMultivariateHistograms(x, ...)
	}

}


# # =================================================================
# # Plot a list of uni.hmms or a matrix with uni.hmm filenames to pdf
# # =================================================================
# plot.distributions.to.pdf <- function(uni.hmms, file='distribution-plots') {
# 
# 	## Preprocess input
# 	uni.hmms <- as.matrix(uni.hmms) # TODO: include handling when only one hmm is given
# 
# 	## Set up the page
# 	library(ggplot2)
# 	library(grid)
# 	ncols <- ncol(uni.hmms)
# 	nrows <- nrow(uni.hmms)
# 	pdf(file=file, width=8*ncols, height=7*nrows)
# 	numPlots <- length(uni.hmms)
# 	grid.newpage()
# 	layout <- matrix(seq(1, ncols * nrows), ncol = ncols, nrow = nrows, byrow = TRUE)
# 	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
# 
# 	## Plot
# 	ptm = proc.time()
# 	for (irow in 1:nrows) {
# 		for (icol in 1:ncols) {
# 			uni.hmm <- uni.hmms[irow,icol]
# 			if (class(uni.hmm)=='character') {
# 				if (uni.hmm!='' & uni.hmm!='NA') {
# 					uni.hmm <- loadHmmsFromFiles(list(uni.hmm))[[1]]
# 					ggplt <- plot.distribution(uni.hmm)
# 				} else {
# 					ggplt <- ggplot(data=data.frame(x=0:10, y=0:10)) + geom_line(aes_string(x='x',y='y'))
# 				}
# 			} else if (class(uni.hmm)==class.univariate.hmm) {
# 				ggplt <- plot.distribution(uni.hmm)
# 			}
# 				
# 			# Get the i,j matrix positions of the regions that contain this subplot
# 			i <- (which(irow==1:nrows)-1) * ncols + which(icol==1:ncols)
# 			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
# 			print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
# 		}
# 	}
# 
# 	d <- dev.off()
# 	print("Time for plotting to file")
# 	print(proc.time() - ptm)
# 
# }

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
plotUnivariateHistogram <- function(model, state=NULL, chrom=NULL, start=NULL, end=NULL) {

	## Check user input
	if (check.univariate.model(model)!=0) {
		model <- get(load(model))
		if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate HMM or a file that contains a univariate HMM")
	}

	## Load libraries
	library(ggplot2)
	library(reshape2)

	# -----------------------------------------
	# Get right x limit
	get_rightxlim <- function(histdata, reads) {
		rightxlim1 <- median(reads[reads>0])*7
		breaks <- histdata$breaks[1:length(histdata$counts)]
		counts <- histdata$counts
		rightxlim2 <- breaks[counts<=5 & breaks>median(reads)][1]
		rightxlim <- min(c(rightxlim1,rightxlim2), na.rm=TRUE)
		return(rightxlim)
	}

	# Select the rows to plot
	selectmask <- rep(TRUE,length(model$bins))
	if (!is.null(chrom)) {
		if (! chrom %in% levels(seqnames(model$bins))) {
			stop(chrom," can't be found in the model coordinates.")
		}
		selectchrom <- as.logical(seqnames(model$bins) == chrom)
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
# Plot a read histogram in normal space of the given state
# ============================================================
plotUnivariateNormalTransformation <- function(model, state='unmodified') {

	## Load libraries
	library(ggplot2)

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

	## Load libraries
	library(ggplot2)

	## Boxplot
	df <- as.data.frame(model$bins)[,c('state','reads')]
	ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes_string(x='state', y='reads', fill='state')) + scale_fill_manual(values=state.colors)
	return(ggplt)

}

# ============================================================
# Plot a heat map of the transition probabilities
# ============================================================
plotMultivariateTransition <- function(multi.hmm) {

	library(ggplot2)
	library(reshape2)

	A <- melt(multi.hmm$transitionProbs, varnames=c('from','to'), value.name='prob')
	A$from <- factor(A$from, levels=stateorderByTransition(multi.hmm))
	A$to <- factor(A$to, levels=stateorderByTransition(multi.hmm))
	ggplt <- ggplot(data=A) + geom_tile(aes_string(x='to', y='from', fill='prob')) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low="white", high="blue")

	return(ggplt)

}
