#' @import ggplot2
NULL

#' Get state colors
#'
#' Get the colors that are used for plotting.
#'
#' @param labels Any combination of \code{c("zero-inflation","unmodified","modified", "total", "counts")}.
#' @return A character vector with colors.
#' @export
#'@examples
#'cols <- getStateColors()
#'pie(1:length(cols), col=cols, labels=names(cols))
getStateColors <- function(labels=NULL) {
    if (is.null(labels)) {
        labels <- c("zero-inflation","unmodified","modified", "total", "counts")
    }
    state.colors <- c("zero-inflation"="gray30","unmodified"="gray48","modified"="orangered3", "total"="black", "counts"="grey35")
    return(state.colors[labels])
}
 
# =================================================================
# Define plotting methods for the generic
# =================================================================
#' Plotting function for saved \pkg{\link{chromstaR}} objects
#'
#' Convenience function that loads and plots a \pkg{\link{chromstaR}} object in one step.
#'
#' @param x A filename that contains either \code{\link{binned.data}}, a \code{\link{uniHMM}} or a \code{\link{multiHMM}}.
#' @param ... Additional arguments.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot character
#' @importFrom graphics plot
#' @export
plot.character <- function(x, ...) {
    x <- get(load(x))
    graphics::plot(x, ...)
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

#' Plotting function for \code{\link{uniHMM}} objects
#'
#' Make different types of plots for \code{\link{uniHMM}} objects.
#'
#' @param x A \code{\link{uniHMM}} object.
#' @param type Type of the plot, one of \code{c('histogram', 'karyogram', 'boxplot', 'normalTransformation')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{histogram}}{A histogram of binned read counts with fitted mixture distribution.}
#'   \item{\code{karyogram}}{A karyogram with binned read counts and peak calls. This uses the \pkg{\link[ggbio]{ggbio}} package and is very slow!}
#'   \item{\code{boxplot}}{A boxplot of read counts for the different states.}
#'   \item{\code{normalTransformation}}{A histogram of transformed read counts.}
#' }
#' @param ... Additional arguments for the different plot types.
#' \describe{
#'   \item{\code{state}}{Plot the \code{histogram}, \code{boxplot} or \code{normalTransformation} only for the specified state. One of \code{c('unmodified','modified')}.}
#'   \item{\code{chromosomes,start,end}}{Plot the \code{histogram} only for the specified chromosomes, start and end position.}
#' }
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot uniHMM
#' @export
plot.uniHMM <- function(x, type='histogram', ...) {
    
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

#' Plotting function for \code{\link{multiHMM}} objects
#'
#' Make different types of plots for \code{\link{multiHMM}} objects.
#'
#' @param x A \code{\link{multiHMM}} object.
#' @param type Type of the plot, one of \code{c('transitionMatrix','histograms','correlation')}. You can also specify the type with an integer number.
#' \describe{
#'   \item{\code{transitionMatrix}}{A heatmap with entries of the transition matrix.}
#'   \item{\code{histograms}}{Fitted histograms of all underlying univariate distributions.}
#'   \item{\code{correlation}}{Correlation between read counts of the different experiments.}
#' }
#' @param ... Additional arguments for the different plot types.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @method plot multiHMM
#' @export
plot.multiHMM <- function(x, type='transitionMatrix', ...) {

    if (type == 'transitionMatrix' | type==1) {
        plotMultivariateTransition(x, ...)
    } else if (type == 'histograms' | type==2) {
        plotMultivariateHistograms(x, ...)
    } else if (type == 'correlation' | type==3) {
        plotMultivariateCorrelation(x, ...)
    } else if (type == 'enrichment' | type==4) {
        plotEnrichment(x, ...)
    }

}

# ============================================================
# Helper functions
# ============================================================
get_rightxlim <- function(counts) {
    if (!any(!counts==0)) {
        return(rightxlim=10)
    }
    rightxlim1 <- median(counts[counts>0])*7
    tab <- table(counts)
    tab <- tab[names(tab)!='0']
    breaks <- as.numeric(names(tab))
    rightxlim2 <- breaks[tab<=5 & breaks>median(counts)*2][1]
    rightxlim <- min(rightxlim1,rightxlim2, na.rm=TRUE)
    if (length(rightxlim)==0 | is.na(rightxlim) | is.infinite(rightxlim)) {
        rightxlim <- 1
    }
    return(rightxlim)
}

# ============================================================
# Plot a read histogram
# ============================================================
#' Plot a histogram of binned read counts
#'
#' Plot a histogram of binned read counts from \code{\link{binned.data}}
#'
#' @param binned.data A \code{\link{binned.data}} object containing binned read counts in meta-column 'counts'.
#' @param chromosomes,start,end Plot the histogram only for the specified chromosomes, start and end position.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom graphics hist
plotBinnedDataHistogram <- function(binned.data, chromosomes=NULL, start=NULL, end=NULL) {

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
    if (length(which(selectmask)) != length(binned.data$counts)) {
        counts <- binned.data$counts[selectmask]
    } else {
        counts <- binned.data$counts
    }


    # Plot the histogram
    rightxlim <- get_rightxlim(counts)
    ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")
    return(ggplt)

}

# =================================================================
# Plot a read histogram with univariate fits for a multivariate HMM
# =================================================================
plotMultivariateHistograms <- function(multi.hmm) {

    ## Make fake uni.hmm and plot
    binmapping <- dec2bin(names(multi.hmm$mapping))
    ggplts <- list()
    for (i1 in 1:length(multi.hmm$IDs)) {
        uni.hmm <- list()
        uni.hmm$ID <- multi.hmm$IDs[i1]
        uni.hmm$bins <- multi.hmm$bins
        mapping <- c('unmodified','modified')[binmapping[,i1]+1]
        names(mapping) <- rownames(binmapping)
        uni.hmm$bins$state <- mapping[uni.hmm$bins$state]
        uni.hmm$bins$counts <- multi.hmm$bins$counts[,i1]
        uni.hmm$weights <- multi.hmm$weights.univariate[[i1]]
        uni.hmm$distributions <- multi.hmm$distributions[[i1]]
        class(uni.hmm) <- class.univariate.hmm
        ggplts[[i1]] <- plotUnivariateHistogram(uni.hmm)
    }
    
    return(ggplts)

}

# =================================================================
# Plot count correlation heatmap for a multivariate HMM
# =================================================================
#' @importFrom stats dist hclust
#' @importFrom reshape2 melt
plotMultivariateCorrelation <- function(multi.hmm) {

    cr <- cor(multi.hmm$bins$counts)
    hc <- stats::hclust(stats::dist(cr))
    df <- reshape2::melt(cr, value.name='correlation')
    df$Var1 <- factor(df$Var1, levels=levels(df$Var1)[hc$order])
    df$Var2 <- factor(df$Var2, levels=levels(df$Var2)[hc$order])
    ggplt <- ggplot(df) + geom_tile(aes_string(x='Var1', y='Var2', fill='correlation')) + xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1)) + scale_fill_gradient(low='red', high='yellow')
    return(ggplt)

}

# ============================================================
# Plot a read histogram with univariate fits
# ============================================================
#' Plot a histogram of binned read counts with fitted mixture distribution
#'
#' Plot a histogram of binned read counts from with fitted mixture distributions from a \code{\link{uniHMM}} object.
#'
#' @param model A \code{\link{uniHMM}} object.
#' @param state Plot the histogram only for the specified state. One of \code{c('unmodified','modified')}.
#' @param chromosomes,start,end Plot the histogram only for the specified chromosomes, start and end position.
#' @param linewidth Width of the distribution lines.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom graphics hist
#' @importFrom stats dnbinom
#' @importFrom reshape2 melt
plotUnivariateHistogram <- function(model, state=NULL, chromosomes=NULL, start=NULL, end=NULL, linewidth=1) {

    ## Check user input
    if (check.univariate.model(model)!=0) {
        model <- get(load(model))
        if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate HMM or a file that contains a univariate HMM")
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
#     if (length(which(selectmask)) != length(model$bins$counts)) {
        counts <- model$bins$counts[selectmask]
        states <- model$bins$state[selectmask]
        weights <- rep(NA, 3)
        weights[1] <- length(which(states=="zero-inflation"))
        weights[2] <- length(which(states=="unmodified"))
        weights[3] <- length(which(states=="modified"))
        weights <- weights / length(states)
#     } else {
#         counts <- model$bins$counts
#         weights <- model$weights
#     }


    # Plot the histogram
    rightxlim <- get_rightxlim(counts)
    ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")

    ### Add fits to the histogram
    x <- 0:max(counts)
    distributions <- data.frame(x)

    # Unmodified
    w <- weights[1]/(weights[2]+weights[1]) # weight for the zero-inflation
    if (is.nan(w)) { w <- 0 }
    distributions$unmodified <- (1-weights[3]) * dzinbinom(x, w, model$distributions[2,'size'], model$distributions[2,'prob'])
    # Modified
    distributions$modified <- weights[3] * stats::dnbinom(x, model$distributions[3,'size'], model$distributions[3,'prob'])
    # Total
    distributions$total <- distributions$unmodified + distributions$modified
    # Convert to long format
    df <- reshape2::melt(distributions, id.vars='x', variable.name='state', value.name='y')

    # Make legend
    lmeans <- round(model$distributions[,'mu'], 2)[-1]
    lvars <- round(model$distributions[,'variance'], 2)[-1]
    lweights <- round(c(1-weights[3], weights[3]), 2)
    legend <- paste0(c('unmodified','modified'), ", mean=", lmeans, ", var=", lvars, ", weight=", lweights)
    legend <- c(legend, paste0('total, mean(data)=', round(mean(counts),2), ', var(data)=', round(var(counts),2)))
    ggplt <- ggplt + ggtitle(model$ID)

    ### Plot the distributions
    if (is.null(state)) {
        ggplt <- ggplt + geom_line(data=df, aes_string(x='x', y='y', col='state'), size=linewidth)
        ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c('unmodified','modified','total')), labels=legend) + theme(legend.justification=c(1,1), legend.position=c(1,1))
    } else {
        if (state=="unmodified") {
            ggplt <- ggplt + geom_line(data=df[df$state=='unmodified',], aes_string(x='x', y='y', col='state'), size=linewidth)
            ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c('unmodified')), labels=legend[1]) + theme(legend.justification=c(1,1), legend.position=c(1,1))
        }
        if (state=="modified") {
            ggplt <- ggplt + geom_line(data=df[df$state=='modified',], aes_string(x='x', y='y', col='state'), size=linewidth)
            ggplt <- ggplt + scale_color_manual(name="components", values=getStateColors(c('modified')), labels=legend[2]) + theme(legend.justification=c(1,1), legend.position=c(1,1))
        }
    }
        
    return(ggplt)

}

# ============================================================
# Plot a karyogram with read counts and univariate calls
# ============================================================
#' Plot a karyogram with read counts and univariate peak calls
#'
#' Plot a karyogram with read counts and peak calls from a \code{\link{uniHMM}} object.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} object or file that contains such an object.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotUnivariateKaryogram <- function(model) {

    ## Check user input
    if (check.univariate.model(model)!=0) {
        model <- get(load(model))
        if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate HMM or a file that contains a univariate HMM")
    }

    # Plot the peaks
    ggplt <- ggbio::autoplot(model$segments[model$segments$state=='modified'], layout='karyogram', color=getStateColors('modified'))

    # Plot the read counts
    ggplt <- ggplt + ggbio::layout_karyogram(model$bins, aes_string(x='start', y='counts'), ylim=c(10,40), geom='line', color=getStateColors('counts'))
    ggplt <- ggplt + theme(axis.text.y=element_blank(), panel.background=element_blank())

    return(ggplt)
}

# ============================================================
# Plot a read histogram in normal space of the given state
# ============================================================
#' @importFrom stats pnbinom qnorm dnorm
plotUnivariateNormalTransformation <- function(model, state='unmodified') {

    ## Plot settings
    cols <- getStateColors(c("unmodified","modified"))

    ## Transform the counts
    df <- as.data.frame(model$bins)[,c('counts','state')]
    # Transform to uniform space
    mask <- df$state=='modified'
    df$ucounts[!mask] <- pzinbinom(df$counts[!mask], model$weights[1], model$distributions[2,'size'], model$distributions[2,'prob'])
    df$ucounts[mask] <- stats::pnbinom(df$counts[mask], model$distributions[3,'size'], model$distributions[3,'prob'])
    # Transform to normal space
    df$ncounts <- stats::qnorm(df$ucounts)

    ## Make the plots
    subset <- df$ncounts[df$state==state]
    breaks <- c(-Inf,sort(as.numeric(names(table(subset)))))
    x <- seq(-4,4,0.1)
    title <- paste0("Transformed emission density for state ",state)
    ggplt <- ggplot() + geom_histogram(data=data.frame(ucounts=subset), aes_string(x='ucounts', y='..density..'), breaks=breaks, right=TRUE, col='black', fill=cols[state]) + theme_bw() + geom_line(data=data.frame(x=x, y=stats::dnorm(x, mean=0, sd=1)), aes_string(x='x', y='y')) + xlab("transformed counts") + labs(title=title)
    return(ggplt)

}

# ============================================================
# Plot a boxplot of the univariate calls
# ============================================================
plotUnivariateBoxplot <- function(model) {

    ## Boxplot
    df <- as.data.frame(model$bins)[,c('state','counts')]
    ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes_string(x='state', y='counts', fill='state')) + scale_fill_manual(values=getStateColors(c('zero-inflation','unmodified','modified')))
    return(ggplt)

}

# ============================================================
# Plot a heat map of the transition probabilities
# ============================================================
#' @importFrom reshape2 melt
plotMultivariateTransition <- function(multi.hmm) {

    A <- reshape2::melt(multi.hmm$transitionProbs, varnames=c('from','to'), value.name='prob')
    A$from <- factor(A$from, levels=stateorderByTransition(multi.hmm))
    A$to <- factor(A$to, levels=stateorderByTransition(multi.hmm))
    ggplt <- ggplot(data=A) + geom_tile(aes_string(x='to', y='from', fill='prob')) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low="white", high="blue")

    return(ggplt)

}

# ============================================================
# Plot a heatmap of combinatorial states
# ============================================================
#' Plot a heatmap of combinatorial states
#'
#' Plot a heatmap that shows the binary presence/absence of marks for the different combinations.
#'
#' @param hmm A \code{\link{multiHMM}} object.
#' @param marks A character vector with histone marks. If specified, \code{hmm} will be ignored.
#' @return A \code{\link[ggplot2]{ggplot}} object containing the plot.
#' @importFrom reshape2 melt
#' @export
#' @author Aaron Taudt
#' @examples 
#'marks <- c('H3K4me3','H3K27me3','H4K20me1')
#'heatmapCombinations(marks=marks)
#'
#'file <- system.file("data","multivariate_mode-mark_condition-SHR.RData",
#'                     package="chromstaR")
#'heatmapCombinations(file)
#'
heatmapCombinations <- function(hmm=NULL, marks=NULL) {

    if (is.null(marks)) {
        hmm <- loadHmmsFromFiles(hmm, check.class=class.multivariate.hmm)[[1]]
        levels.combinations <- levels(hmm$bins$combination)
        levels.combinations <- gsub('\\[', '', levels.combinations)
        levels.combinations <- gsub('\\]', '', levels.combinations)
        marks <- unique(unlist(strsplit(levels.combinations, '\\+')))
    }
    d <- dec2bin(0:(2^length(marks)-1), colnames=marks)
    d <- as.data.frame(d)
    d$combination <- apply(d, 1, function(x) { paste(colnames(d)[x],collapse='+') })
    d$combination <- paste0('[',d$combination,']')
    df <- reshape2::melt(d, variable.name='mark', value.name='emission', id.vars='combination')
    df$emission <- as.numeric(df$emission)

    ggplt <- ggplot(df) + geom_tile(aes_string(x='combination', y='mark', fill='emission')) + scale_fill_gradient(low='white', high='blue', guide=FALSE, limits=c(0,1)) + theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))

    return(ggplt)
}

