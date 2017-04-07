#' @import ggplot2
NULL

#' Get state colors
#'
#' Get the colors that are used for plotting.
#'
#' @param labels Any combination of \code{c("zero-inflation","unmodified","modified", "total", "counts")}.
#' @return A character vector with colors.
#' @seealso \code{\link{plotting}}
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
# Multivariate plotting functions
# ============================================================
#' Histograms of binned read counts with fitted mixture distribution
#'
#' Plot histograms of binned read counts with fitted mixture distributions from a \code{\link{multiHMM}} object.
#'
#' @param model A \code{\link{multiHMM}} object or file that contains such an object.
#' @param ... Additional arguments (see \code{\link{plotHistogram}}).
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @seealso \code{\link{plotting}}
plotHistograms <- function(model, ...) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=class.multivariate.hmm)[[1]] )
    ## Make fake uni.hmm and plot
    binmapping <- dec2bin(levels(model$bins$state), colnames=model$info$ID)
    ggplts <- list()
    for (i1 in 1:ncol(model$bins$counts)) {
        uni.hmm <- list()
        uni.hmm$info <- model$info[i1,]
        uni.hmm$bins <- model$bins
        mapping <- c('unmodified','modified')[binmapping[,i1]+1]
        names(mapping) <- rownames(binmapping)
        uni.hmm$bins$state <- mapping[uni.hmm$bins$state]
        uni.hmm$bins$counts <- model$bins$counts[,i1]
        uni.hmm$weights <- model$weights.univariate[[i1]]
        uni.hmm$distributions <- model$distributions[[i1]]
        class(uni.hmm) <- class.univariate.hmm
        ggplts[[i1]] <- plotHistogram(uni.hmm, ...)
    }
    
    return(ggplts)

}

#' Read count correlation heatmap
#'
#' Heatmap of read count correlations (see \code{\link[stats]{cor}}).
#' 
#' @param model A \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object or file that contains such an object.
#' @param cluster Logical indicating whether or not to cluster the heatmap.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @importFrom stats dist hclust
#' @importFrom reshape2 melt
#' @seealso \code{\link{plotting}}
#' @export
#' @examples 
#'## Get an example multiHMM ##
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'## Plot count correlations as heatmap
#'heatmapCountCorrelation(model)
#'
heatmapCountCorrelation <- function(model, cluster=TRUE) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]] )
    cr <- cor(model$bins$counts.rpkm)
    df <- reshape2::melt(cr, value.name='correlation')
    if (cluster) {
        hc <- stats::hclust(stats::dist(cr))
        df$Var1 <- factor(df$Var1, levels=levels(df$Var1)[hc$order])
        df$Var2 <- factor(df$Var2, levels=levels(df$Var2)[hc$order])
    }
    ggplt <- ggplot(df) + geom_tile(aes_string(x='Var1', y='Var2', fill='correlation')) + xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + scale_fill_gradient(low='red', high='yellow')
    return(ggplt)

}

#' Heatmap of transition probabilities
#'
#' Plot a heatmap of transition probabilities for a \code{\link{multiHMM}} model.
#'
#' @param model A \code{\link{multiHMM}} object or file that contains such an object.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @importFrom reshape2 melt
#' @seealso \code{\link{plotting}}
#' @export
#' @examples 
#'## Get an example multiHMM ##
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'## Plot transition probabilites as heatmap
#'heatmapTransitionProbs(model)
#'
heatmapTransitionProbs <- function(model) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=class.multivariate.hmm)[[1]] )
    A <- reshape2::melt(model$transitionProbs, varnames=c('from','to'), value.name='prob')
    A$from <- factor(A$from, levels=stateorderByTransition(model))
    A$to <- factor(A$to, levels=stateorderByTransition(model))
    ggplt <- ggplot(data=A) + geom_tile(aes_string(x='to', y='from', fill='prob')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_gradient(low="white", high="blue")

    return(ggplt)

}

#' Plot a heatmap of combinatorial states
#'
#' Plot a heatmap that shows the binary presence/absence of marks for the different combinations.
#'
#' @param model A \code{\link{multiHMM}} object or file that contains such an object.
#' @param marks A character vector with histone marks. If specified, \code{model} will be ignored.
#' @param emissionProbs A matrix with emission probabilities where \code{dimnames(emissionProbs)} gives the state labels and marks. This option is helpful to plot probabilistic chromatin states (not part of \pkg{\link{chromstaR}}). If specified, \code{model} and \code{marks} will be ignored.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @importFrom reshape2 melt
#' @seealso \code{\link{plotting}}
#' @export
#' @author Aaron Taudt
#' @examples 
#'marks <- c('H3K4me3','H3K27me3','H4K20me1')
#'heatmapCombinations(marks=marks)
#'
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'heatmapCombinations(file)
#'
heatmapCombinations <- function(model=NULL, marks=NULL, emissionProbs=NULL) {

    if (!is.null(emissionProbs)) {
        df <- reshape2::melt(emissionProbs, value.name='emission')
        names(df)[1:2] <- c('combination', 'mark')
        df$combination <- factor(df$combination)
    } else {
        if (is.null(marks)) {
            model <- loadHmmsFromFiles(model, check.class=class.multivariate.hmm)[[1]]
            levels.combinations <- levels(model$bins$combination)
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
    }

    ggplt <- ggplot(df) + geom_tile(aes_string(x='combination', y='mark', fill='emission')) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
    if (!is.null(emissionProbs)) {
        ggplt <- ggplt + scale_fill_gradient(low='white', high='blue', limits=c(0,1))
    } else {
        ggplt <- ggplt + scale_fill_gradient(low='white', high='blue', guide=FALSE, limits=c(0,1))
    }

    return(ggplt)
}

# ============================================================
# Univariate plotting functions
# ============================================================
#' Histogram of binned read counts with fitted mixture distribution
#'
#' Plot a histogram of binned read counts with fitted mixture distributions from a \code{\link{uniHMM}} object.
#'
#' @param model A \code{\link{uniHMM}} object or file that contains such an object.
#' @param state Plot the histogram only for the specified state. One of \code{c('unmodified','modified')}.
#' @param chromosomes,start,end Plot the histogram only for the specified chromosomes, start and end position.
#' @param linewidth Width of the distribution lines.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom graphics hist
#' @importFrom stats dnbinom
#' @importFrom reshape2 melt
#' @seealso \code{\link{plotting}}
#' @export
#' @examples
#'## Get an example BAM file with ChIP-seq reads
#'file <- system.file("extdata", "euratrans",
#'                       "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                        package="chromstaRData")
#'## Bin the BED file into bin size 1000bp
#'data(rn4_chrominfo)
#'data(experiment_table)
#'binned <- binReads(file, experiment.table=experiment_table,
#'                   assembly=rn4_chrominfo, binsizes=1000,
#'                   stepsizes=500, chromosomes='chr12')
#'plotHistogram(binned)
#'## Fit the univariate Hidden Markov Model
#'hmm <- callPeaksUnivariate(binned, max.time=60, eps=1)
#'## Check if the fit is ok
#'plotHistogram(hmm)
#'
plotHistogram <- function(model, state=NULL, chromosomes=NULL, start=NULL, end=NULL, linewidth=1) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=c('GRanges', class.univariate.hmm))[[1]] )
    if (is(model, 'GRanges')) {
        bins <- model
    } else if (is(model, class.univariate.hmm)) {
        bins <- model$bincounts
    }

    # Select the rows to plot
    selectmask <- rep(TRUE,length(bins))
    if (!is.null(chromosomes)) {
        if (any(! chromosomes %in% levels(seqnames(bins)))) {
            stop(chromosomes[! chromosomes %in% levels(seqnames(bins))]," can't be found in the binned data.")
        }
        selectchrom <- as.logical(seqnames(bins) %in% chromosomes)
        selectmask <- selectmask & selectchrom
        if (!is.null(start)) {
            selectstart <- as.logical(start(ranges(bins)) >= start)
            selectmask <- selectmask & selectstart
        }
        if (!is.null(end)) {
            selectend <- as.logical(end(ranges(bins)) <= end)
            selectmask <- selectmask & selectend
        }
    }
    if (!is.null(state)) {
        selectmask <- selectmask & bins$state==state
    }
    counts <- bins$counts[selectmask, '0']

    # Plot the histogram
    rightxlim <- get_rightxlim(counts)
    ggplt <- ggplot(data.frame(counts)) + geom_histogram(aes_string(x='counts', y='..density..'), binwidth=1, color='black', fill='white') + coord_cartesian(xlim=c(0,rightxlim)) + theme_bw() + xlab("read count")
    
    if (is(model, 'GRanges')) {
        return(ggplt)
    }

    ### Add fits to the histogram
    states <- model$bins$state[selectmask]
    weights <- rep(NA, 3)
    weights[1] <- length(which(states=="zero-inflation"))
    weights[2] <- length(which(states=="unmodified"))
    weights[3] <- length(which(states=="modified"))
    weights <- weights / length(states)

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
    ggplt <- ggplt + ggtitle(model$info$ID)

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

#' Plot a karyogram with read counts and univariate peak calls
#'
#' Plot a karyogram with read counts and peak calls from a \code{\link{uniHMM}} object.
#'
#' @author Aaron Taudt
#' @param model A \code{\link{uniHMM}} object or file that contains such an object.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotKaryogram <- function(model) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=class.univariate.hmm)[[1]] )

    # Plot the peaks
    ggplt <- ggbio::autoplot(model$segments[model$segments$state=='modified'], layout='karyogram', color=getStateColors('modified'))

    # Plot the read counts
    ggplt <- ggplt + ggbio::layout_karyogram(model$bins, aes_string(x='start', y='counts'), ylim=c(10,40), geom='line', color=getStateColors('counts'))
    ggplt <- ggplt + theme(axis.text.y=element_blank(), panel.background=element_blank())

    return(ggplt)
}

#' @importFrom stats pnbinom qnorm dnorm
plotNormalTransformation <- function(model, state='unmodified') {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=class.univariate.hmm)[[1]] )

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

plotBoxplot <- function(model) {

    model <- suppressMessages( loadHmmsFromFiles(model, check.class=class.univariate.hmm)[[1]] )

    ## Boxplot
    df <- as.data.frame(model$bins)[,c('state','counts')]
    ggplt <- ggplot() + theme_bw() + geom_boxplot(data=df, aes_string(x='state', y='counts', fill='state')) + scale_fill_manual(values=getStateColors(c('zero-inflation','unmodified','modified')))
    return(ggplt)

}


#' Plot a genome browser view
#' 
#' Plot a simple genome browser view. This is useful for scripted genome browser snapshots.
#' 
#' @param counts A \code{\link[GenomicRanges]{GRanges}} object with meta-data column 'counts'.
#' @param peaklist A named list() of \code{\link[GenomicRanges]{GRanges}} objects containing peak coordinates.
#' @param chr,start,end Chromosome, start and end coordinates for the plot.
#' @param countcol A character giving the color for the counts.
#' @param peakcols A character vector with colors for the peaks in \code{peaklist}.
#' @param style One of \code{c('peaks', 'density')}.
#' @param peakTrackHeight Relative height of the tracks given in \code{peaklist} compared to the \code{counts}.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
plotGenomeBrowser <- function(counts, peaklist=NULL, chr, start, end, countcol='black', peakcols=NULL, style='peaks', peakTrackHeight=5) {
  
    ## Select ranges to plot
    ranges2plot <- reduce(counts[counts@seqnames == chr & start(counts) >= start & start(counts) <= end])
    
    ## Counts
    counts <- subsetByOverlaps(counts, ranges2plot)
    
    if (style == 'peaks') {
        # df.start <- data.frame(x=start(counts), counts=counts$counts)
        # df.end <- data.frame(x=end(counts), counts=counts$counts)
        # df <- rbind(df.start, df.end)
        # df <- df[order(df$x),]
        df <- data.frame(x=(start(counts)+end(counts))/2, counts=counts$counts) # plot triangles centered at middle of the bin
        ggplt <- ggplot(df) + geom_area(aes_string(x='x', y='counts')) + theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank())
        maxcounts <- max(counts$counts)
        ggplt <- ggplt + scale_y_continuous(breaks=c(0, maxcounts))
    } else if (style == 'density') {
        df <- data.frame(xmin=start(counts), xmax=end(counts), counts=counts$counts)
        
        # # Rolling mean
        # n <- 10
        # cx <- cumsum(df$counts)
        # rsum <- (cx[n:length(df$counts)] - c(0, cx[1:(length(df$counts) - n)])) / n
        # df$counts[(n%/%2 + 1):(length(df$counts)-(n-1)%/%2)] <- rsum
        
        # # Expand high peaks
        # fact <- 1e-5
        # df$xmin <- df$xmin - (end-start)*df$counts * fact
        # df$xmax <- df$xmax + (end-start)*df$counts * fact
        
        ggplt <- ggplot(df) + geom_rect(aes_string(xmin='xmin', xmax='xmax', ymin=0, ymax=4, alpha='counts')) + theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    } else {
        stop("Unknown value '", style, "' for parameter 'style'. Must be one of c('peaks', 'density').")
    }
    
    ## Peaks
    if (!is.null(peaklist)) {
        if (is.null(peakcols)) {
            peakcols <- getDistinctColors(length(peaklist))
        }
        for (i1 in 1:length(peaklist)) {
            p <- peakTrackHeight
            peaks <- subsetByOverlaps(peaklist[[i1]], ranges2plot)
            if (length(peaks) > 0) {
                df <- data.frame(start=start(peaks), end=end(peaks), ymin=-p*i1, ymax=-p*i1+0.9*p)
                ggplt <- ggplt + geom_rect(data=df, mapping=aes_string(xmin='start', xmax='end', ymin='ymin', ymax='ymax'), col=peakcols[i1], fill=peakcols[i1])
            }
            trackname <- names(peaklist)[i1]
            df <- data.frame(x=start(counts)[1], y=-p*i1+0.5*p, label=trackname)
            ggplt <- ggplt + geom_text(data=df, mapping=aes_string(x='x', y='y', label='label'), vjust=0.5, hjust=0.5, col=peakcols[i1])
        }
    }
    
    return(ggplt)
}