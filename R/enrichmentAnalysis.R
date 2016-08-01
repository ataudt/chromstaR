#' Enrichment analysis
#' 
#' Plotting functions for enrichment analysis of \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} objects with any annotation of interest, specified as a \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @name enrichment_analysis
#' @param combinations A vector with combinations for which the fold enrichment will be calculated. If \code{NULL} all combinations will be considered.
#' @param marks A vector with marks for which the enrichment is plotted. If \code{NULL} all marks will be considered.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object containing the plot or a list() with \code{\link[ggplot2:ggplot]{ggplot}} objects if several plots are returned. For \code{plotFoldEnrichHeatmap} a named array with fold enrichments if \code{plot=FALSE}.
#' @author Aaron Taudt
#' @seealso \code{\link{plotting}}
#' @examples 
#'### Get an example multiHMM ###
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'
#'### Obtain gene coordinates for rat from biomaRt ###
#'library(biomaRt)
#'ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host='may2012.archive.ensembl.org',
#'                   dataset='rnorvegicus_gene_ensembl')
#'genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position',
#'                            'end_position', 'strand', 'external_gene_id',
#'                            'gene_biotype'),
#'               mart=ensembl)
#'# Transform to GRanges for easier handling
#'genes <- GRanges(seqnames=paste0('chr',genes$chromosome_name),
#'                 ranges=IRanges(start=genes$start, end=genes$end),
#'                 strand=genes$strand,
#'                 name=genes$external_gene_id, biotype=genes$gene_biotype)
#'print(genes)
#'
#'### Make the enrichment plots ###
#'# We expect promoter [H3K4me3] and bivalent-promoter signatures [H3K4me3+H3K27me3]
#'# to be enriched at transcription start sites.
#'    plotFoldEnrichment(hmm = model, annotation = genes, bp.around.annotation = 15000) +
#'    ggtitle('Fold enrichment around genes') +
#'    xlab('distance from gene body')
#'  
#'# Plot enrichment only at TSS. We make use of the fact that TSS is the start of a gene.
#'    plotFoldEnrichment(model, genes, region = 'start') +
#'    ggtitle('Fold enrichment around TSS') +
#'    xlab('distance from TSS in [bp]')
#'# Note: If you want to facet the plot because you have many combinatorial states you
#'# can do that with
#'    plotFoldEnrichment(model, genes, region = 'start') +
#'    facet_wrap(~ combination)
#'  
#'# Another form of visualization that shows every TSS in a heatmap
#'# If transparency is not supported try to plot to pdf() instead.
#'    tss <- resize(genes, width = 3, fix = 'start')
#'    plotEnrichCountHeatmap(model, tss) +
#'    theme(strip.text.x = element_text(size=6))
#'  
#'# Fold enrichment with different biotypes, showing that protein coding genes are
#'# enriched with (bivalent) promoter combinations [H3K4me3] and [H3K4me3+H3K27me3],
#'# while rRNA is enriched with the empty [] and repressive combinations [H3K27me3].
#'    biotypes <- split(tss, tss$biotype)
#'    plotFoldEnrichHeatmap(model, annotations=biotypes) + coord_flip()
#'
NULL


#' @describeIn enrichment_analysis Compute the fold enrichment of combinatorial states for multiple annotations.
#' @param hmm A \code{\link{combinedMultiHMM}} or \code{\link{multiHMM}} object or a file that contains such an object.
#' @param annotations A \code{list()} with \code{\link{GRanges}} objects containing coordinates of multiple annotations The names of the list entries will be used to name the return values.
#' @param plot A logical indicating whether the plot or an array with the fold enrichment values is returned.
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom reshape2 melt
#' @export
plotFoldEnrichHeatmap <- function(hmm, annotations, what="combinations", combinations=NULL, marks=NULL, plot=TRUE) {
    
    hmm <- loadHmmsFromFiles(hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    ## Variables
    bins <- hmm$bins
    if (class(hmm) == class.combined.multivariate.hmm) {
    } else if (class(hmm) == class.multivariate.hmm) {
        # Rename 'combination' to 'combination.' for coherence with combinedMultiHMM
        names(mcols(bins))[grep('combination', names(mcols(bins)))] <- 'combination.'
    }
    conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
    if (is.null(combinations)) {
        comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
    } else {
        comb.levels <- combinations
    }
    if (is.null(marks)) {
        mark.levels <- unique(hmm$info$mark)
    } else {
        mark.levels <- marks
    }
    genome <- sum(as.numeric(width(bins)))
    feature.lengths <- lapply(annotations, function(x) { sum(as.numeric(width(x))) })
    
    if (what == 'peaks') {
        binstates <- dec2bin(bins$state, colnames=hmm$info$ID)
    }
    
    ## Fold enrichment
    ggplts <- list()
    folds <- list()
    for (condition in conditions) {
        if (what == 'combinations') {
            bins$combination <- mcols(bins)[,paste0('combination.', condition)]
            fold <- array(NA, dim=c(length(annotations), length(comb.levels)), dimnames=list(annotation=names(annotations), combination=comb.levels))
            for (icomb in 1:length(comb.levels)) {
                mask <- bins$combination == comb.levels[icomb]
                bins.mask <- bins[mask]
                combstate.length <- sum(as.numeric(width(bins.mask)))
                for (ifeat in 1:length(annotations)) {
                    feature <- annotations[[ifeat]]
                    ind <- findOverlaps(bins.mask, feature)
    
                    binsinfeature <- bins.mask[unique(S4Vectors::queryHits(ind))]
                    sum.binsinfeature <- sum(as.numeric(width(binsinfeature)))
    
                    featuresinbins <- feature[unique(S4Vectors::subjectHits(ind))]
                    sum.featuresinbins <- sum(as.numeric(width(featuresinbins)))
    
                    fold[ifeat,icomb] <- sum.binsinfeature / combstate.length / feature.lengths[[ifeat]] * genome
                }
            }
    
        } else if (what == 'peaks') {
            fold <- array(NA, dim=c(length(annotations), length(mark.levels)), dimnames=list(annotation=names(annotations), mark=mark.levels))
            for (imark in 1:length(mark.levels)) {
                mark <- mark.levels[imark]
                colmask <- hmm$info$mark == mark
                colmask <- colmask & (!duplicated(paste0(hmm$info$mark, hmm$info$condition))) # remove replicates
                if (condition != "") {
                    colmask <- colmask & (hmm$info$condition == condition)
                }
                binstates.cond <- binstates[,colmask]
                mask <- binstates.cond
                bins.mask <- bins[mask]
                combstate.length <- sum(as.numeric(width(bins.mask)))
                for (ifeat in 1:length(annotations)) {
                    feature <- annotations[[ifeat]]
                    ind <- findOverlaps(bins.mask, feature)
    
                    binsinfeature <- bins.mask[unique(S4Vectors::queryHits(ind))]
                    sum.binsinfeature <- sum(as.numeric(width(binsinfeature)))
    
                    featuresinbins <- feature[unique(S4Vectors::subjectHits(ind))]
                    sum.featuresinbins <- sum(as.numeric(width(featuresinbins)))
    
                    fold[ifeat,imark] <- sum.binsinfeature / combstate.length / feature.lengths[[ifeat]] * genome
                }
            }
          
        }
      
        if (plot) {
            df <- reshape2::melt(fold, value.name='foldEnrichment')
            if (what == 'combinations') {
                ggplt <- ggplot(df) + geom_tile(aes_string(x='combination', y='annotation', fill='foldEnrichment')) + theme(axis.text.x = element_text(angle=90, hjust=1)) + scale_fill_gradient(low='white', high='blue')
            } else if (what == 'peaks') {
                ggplt <- ggplot(df) + geom_tile(aes_string(x='mark', y='annotation', fill='foldEnrichment')) + theme(axis.text.x = element_text(angle=90, hjust=1)) + scale_fill_gradient(low='white', high='blue')
            }
            ggplts[[condition]] <- ggplt
        } else {
            folds[[condition]] <- fold
        }
    }

    if (class(hmm) == class.multivariate.hmm) {
        if (plot) {
            return(ggplts[[1]])
        } else {
            return(folds[[1]])
        }
    } else if (class(hmm) == class.combined.multivariate.hmm) {
        if (plot) {
            return(ggplts)
        } else {
            return(folds)
        }
    }
    
}


#' @describeIn enrichment_analysis Plot read counts around annotation as heatmap.
#' @inheritParams enrichmentAtAnnotation
#' @param max.rows An integer specifying the number of randomly subsampled rows that are plotted from the \code{annotation} object. This is necessary to avoid crashing for heatmaps with too many rows.
#' @importFrom reshape2 melt
#' @export
plotEnrichCountHeatmap <- function(hmm, annotation, bp.around.annotation=10000, max.rows=1000) {

    hmm <- loadHmmsFromFiles(hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    ## Variables
    bins <- hmm$bins
    if (class(hmm) == class.combined.multivariate.hmm) {
    } else if (class(hmm) == class.multivariate.hmm) {
        # Rename 'combination' to 'combination.' for coherence with combinedMultiHMM
        names(mcols(bins))[grep('combination', names(mcols(bins)))] <- 'combination.'
    }
    conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
    comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
    binsize <- width(bins)[1]
    around <- round(bp.around.annotation/binsize)
    ## Create new column combination with all conditions combined
    combinations <- list()
    for (condition in conditions) {
        combinations[[condition]] <- paste0(condition, ":", mcols(bins)[,paste0('combination.', condition)])
    }
    combinations$sep <- ', '
    bins$combination <- factor(do.call(paste, combinations))

    # Subsampling for plotting of huge data.frames
    annotation <- subsetByOverlaps(annotation, bins)
    if (length(annotation)>max.rows) {
        annotation <- sample(annotation, size=max.rows, replace=FALSE)
    }
  
    # Get bins that overlap the start of annotation
    ptm <- startTimedMessage("Overlaps with annotation ...")
    index.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="first")
    index.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="last")
    index.plus <- index.plus[!is.na(index.plus)]
    index.minus <- index.minus[!is.na(index.minus)]
    index <- c(index.plus, index.minus)
    stopTimedMessage(ptm)
    
    # Get surrounding indices
    ptm <- startTimedMessage("Getting surrounding indices ...")
    ext.index.plus <- array(NA, dim=c(length(index.plus), 2*around+1), dimnames=list(anno=1:length(index.plus), position=binsize*seq(-around, around, 1)))
    for (i1 in 1:length(index.plus)) {
        ext.index.plus[i1,] <- seq(from=-around+index.plus[i1], to=index.plus[i1]+around)
    }
    if (length(index.minus)>0) {
        ext.index.minus <- array(NA, dim=c(length(index.minus), 2*around+1), dimnames=list(anno=1:length(index.minus), position=binsize*seq(-around, around, 1)))
        for (i1 in 1:length(index.minus)) {
            ext.index.minus[i1,] <- rev( seq(from=-around+index.minus[i1], to=index.minus[i1]+around) )
        }
        ext.index <- rbind(ext.index.plus, ext.index.minus)
    } else {
        ext.index <- ext.index.plus
    }
    ext.index[ext.index <= 0] <- NA
    ext.index[ext.index > length(bins)] <- NA
    rownames(ext.index) <- 1:nrow(ext.index)
    stopTimedMessage(ptm)
    
    ## Go through combinations and then tracks to get the read counts
    ptm <- startTimedMessage("Getting read counts")
    counts <- list()
    for (combination in levels(bins$combination)) {
        counts[[combination]] <- list()
        index.combination <- which(bins$combination[index]==combination)
        ext.index.combination <- ext.index[index.combination,]
        if (is.null(dim(ext.index.combination))) {
            ext.index.combination <- array(ext.index.combination, dim=c(1,dim(ext.index)[[2]]), dimnames=list(anno=rownames(ext.index)[index.combination], position=dimnames(ext.index)[[2]]))
        }
        for (ntrack in colnames(bins$counts)) {
            counts[[combination]][[ntrack]] <- array(bins$counts[ext.index.combination,ntrack], dim=dim(ext.index.combination), dimnames=dimnames(ext.index.combination))
        }
    }
    stopTimedMessage(ptm)
    
    ## Theme
    custom_theme <- theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill='NA'),
        panel.background = element_rect(fill='white'),
        axis.text.y = element_blank()
    )
    
    ## Prepare data.frame
    ptm <- startTimedMessage("Making the plot ...")
    # Exclude rare combinations for plotting
    num.comb <- sapply(counts, function(x) { nrow(x[[1]]) })
    comb2keep <- names(num.comb)[num.comb/sum(num.comb) > 0.005]
    counts <- counts[comb2keep]
    df <- reshape2::melt(counts)
    names(df) <- c('id','position','counts','track','combination')
    df$id <- factor(df$id, levels=rev(unique(df$id)))
    df$combination <- factor(df$combination, levels=unique(df$combination))
    df$track <- factor(df$track, levels=colnames(bins$counts))
    
    ## Plot as heatmap
    ggplt <- ggplot(df) + geom_tile(aes_string(x='position', y='id', color='combination'))
    ggplt <- ggplt + geom_tile(aes_string(x='position', y='id', fill='counts'), alpha=0.6)
    ggplt <- ggplt + facet_wrap( ~ track, nrow=1) + custom_theme
    ggplt <- ggplt + xlab('distance from annotation in [bp]') + ylab('')
    breaks <- sort(c(0, 10^(0:5), max(df$counts, na.rm = TRUE)))
    ggplt <- ggplt + scale_fill_continuous(trans='log1p', breaks=breaks, labels=breaks, low='white', high='black')
    # Insert horizontal lines
    y.lines <- sapply(split(df$id, df$combination), function(x) { max(as.integer(x)) })
    df.lines <- data.frame(y=sort(y.lines[-1]) + 0.5)
    ggplt <- ggplt + geom_hline(data=df.lines, mapping=aes_string(yintercept='y'), linetype=2)
    # Increase color size in legend
    ggplt <- ggplt + guides(color=guide_legend(override.aes = list(size=2)))
    stopTimedMessage(ptm)
    
    return(ggplt)
}


#' @describeIn enrichment_analysis Plot fold enrichment of combinatorial states around and inside of annotation.
#' @inheritParams enrichmentAtAnnotation
#' @importFrom reshape2 melt
#' @export
plotFoldEnrichment <- function(hmm, annotation, bp.around.annotation=10000, region=c("start","inside","end"), num.intervals=20, what='combinations', combinations=NULL, marks=NULL) {

    ## Check user input
    if ((!what %in% c('combinations','peaks','counts')) | length(what) > 1) {
        stop("argument 'what' must be one of c('combinations','peaks','counts')")
    }
  
    ## Variables
    hmm <- loadHmmsFromFiles(hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    bins <- hmm$bins
    if (class(hmm) == class.combined.multivariate.hmm) {
    } else if (class(hmm) == class.multivariate.hmm) {
        # Rename 'combination' to 'combination.' for coherence with combinedMultiHMM
        names(mcols(bins))[grep('combination', names(mcols(bins)))] <- 'combination.'
    }
    conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
    if (is.null(combinations)) {
        comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
    } else {
        comb.levels <- combinations
    }
    if (is.null(marks)) {
        mark.levels <- unique(hmm$info$mark)
    } else {
        mark.levels <- marks
    }

    if (what %in% c('peaks','counts')) {
        ### Get fold enrichment
        enrich <- enrichmentAtAnnotation(hmm$bins, hmm$info, annotation, bp.around.annotation=bp.around.annotation, region=region, what=what, num.intervals=num.intervals)
    }
    ggplts <- list()
    for (condition in conditions) {
        if (what == 'combinations') {
            ### Get fold enrichment
            bins$combination <- mcols(bins)[,paste0('combination.', condition)]
            enrich.cond <- enrichmentAtAnnotation(bins, hmm$info, annotation, bp.around.annotation=bp.around.annotation, region=region, what=what, num.intervals=num.intervals)
        } else {
            enrich.cond <- enrich
        }
        ### Prepare for plotting
        df <- reshape2::melt(enrich.cond)
        df$L1 <- factor(df$L1, levels=c('start','inside','end'))
        df <- rbind(df[df$L1 == 'start',], df[df$L1 == 'inside',], df[df$L1 == 'end',])
        if (length(region)>=2) {
            df <- df[!(df$L1 == 'start' & df$lag > 0),]
            df <- df[!(df$L1 == 'end' & df$lag < 0),]
            df$position <- apply(data.frame(df$interval, df$lag), 1, max, na.rm = TRUE)
        } else if (length(region)==1) {
            df <- df[df$L1 == region,]
            df$position <- df$lag
        }
        df$position[df$L1 == 'end'] <- df$position[df$L1 == 'end'] + bp.around.annotation
        df$position[df$L1 == 'inside'] <- df$position[df$L1 == 'inside'] * bp.around.annotation
        if (what == 'combinations') {
            df <- df[df$combination %in% comb.levels,]
        } else if (what %in% c('peaks','counts')) {
            df$mark <- sub("-.*", "", df$track)
            df <- df[df$mark %in% mark.levels, ]
            df$condition <- sapply(strsplit(as.character(df$track), '-'), '[[', 2)
            if (condition != "") {
                df <- df[df$condition == condition, ]
            }
        }

        ### Plot
        if (what == 'combinations') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='combination'), size=2) + ylab('fold enrichment')
        } else if (what == 'peaks') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='mark'), size=2) + ylab('fraction of positions in peak')
        } else if (what == 'counts') {
            ggplt <- ggplot(df) + geom_line(aes_string(x='position', y='value', col='track'), size=2) + ylab('read count')
        }
        ggplt <- ggplt + theme_bw() + xlab('distance from annotation in [bp]')
        if (length(region)>=2) {
            breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * bp.around.annotation)
            labels <- c(-bp.around.annotation, -bp.around.annotation/2, '0%', '50%', '100%', bp.around.annotation/2, bp.around.annotation)
            ggplt <- ggplt + scale_x_continuous(breaks=breaks, labels=labels)
        }
        ggplts[[condition]] <- ggplt
    }
    if (class(hmm) == class.multivariate.hmm) {
        return(ggplts[[1]])
    } else if (class(hmm) == class.combined.multivariate.hmm) {
        return(ggplts)
    }
    
}


#' Enrichment of (combinatorial) states for genomic annotations
#'
#' The function calculates the enrichment of a genomic feature with peaks or combinatorial states. Input is a \code{\link{multiHMM}} object (containing the peak calls and combinatorial states) and a \code{\link{GRanges}} object containing the annotation of interest (e.g. transcription start sites or genes).
#'
#' @author Aaron Taudt
#' @param bins The \code{$bins} entry from a \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object.
#' @param info The \code{$info} entry from a \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object.
#' @param annotation A \code{\link{GRanges}} object with the annotation of interest.
#' @param bp.around.annotation An integer specifying the number of basepairs up- and downstream of the annotation for which the enrichment will be calculated.
#' @param region A combination of \code{c('start','inside','end')} specifying the region of the annotation for which the enrichment will be calculated. Select \code{'start'} if you have a point-sized annotation like transcription start sites. Select \code{c('start','inside','end')} if you have long annotations like genes.
#' @param what One of \code{c('combinations','peaks','counts')} specifying which statistic to calculate.
#' @param num.intervals Number of intervals for enrichment 'inside' of annotation.
#' @return A \code{list()} containing \code{data.frame()}s for enrichment of combinatorial states and binary states at the start, end and inside of the annotation.
#' @importFrom S4Vectors as.factor subjectHits queryHits
enrichmentAtAnnotation <- function(bins, info, annotation, bp.around.annotation=10000, region=c('start','inside','end'), what=c('combinations','peaks','counts'), num.intervals=21) {

    ## Check user input
    if ((!what %in% c('combinations','peaks','counts')) | length(what) > 1) {
        stop("argument 'what' must be one of c('combinations','peaks','counts')")
    }
  
    ## Variables
    binsize <- width(bins)[1]
    lag <- round(bp.around.annotation/binsize)
    enrich <- list()
    enrich$combinations <- list()
    enrich$peaks <- list()
    enrich$counts <- list()
    info.dedup <- info[!duplicated(paste0(info$mark, info$condition)), ]

    ## Get combinatorial and binary states
    combinations <- bins$combination
    tcombinations <- table(combinations)
    if ('peaks' %in% what) {
        binstates <- dec2bin(bins$state, colnames=info$ID)
        # Remove replicates
        binstates <- binstates[ ,info.dedup$ID]
    }
    if ('counts' %in% what) {
        counts <- bins$counts
    }

    ### Calculating enrichment inside of annotation ###
    if ('inside' %in% region) {
        ptm <- startTimedMessage("Enrichment inside of annotations ...")

        intervals <- seq(from=0, to=1, length.out=num.intervals+1)
        widths.annotation <- width(annotation) - 1
        annotation.1bp <- resize(annotation, 1, fix='start')
        # Initialize arrays
        if ('peaks' %in% what) binstates.inside <- array(dim=c(num.intervals+1, length(info.dedup$ID)), dimnames=list(interval=intervals, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.inside <- array(dim=c(num.intervals+1, length(levels(bins$combination))), dimnames=list(interval=intervals, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.inside <- array(dim=c(num.intervals+1, length(info$ID)), dimnames=list(interval=intervals, track=info$ID))

        for (interval in intervals) {
            shift <- widths.annotation * interval * c(1,-1,1)[as.integer(strand(annotation))]
            shifted.starts <- start(annotation.1bp) + shift
            annotation.shifted <- GRanges(seqnames = seqnames(annotation.1bp), ranges = IRanges(start = shifted.starts, end = shifted.starts), strand = strand(annotation.1bp))
            # Get bins that overlap the shifted annotation
            index.inside.plus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='+' | strand(annotation.shifted)=='*'], bins, select="first")
            index.inside.minus <- findOverlaps(annotation.shifted[strand(annotation.shifted)=='-'], bins, select="last")
            index.inside.plus <- index.inside.plus[!is.na(index.inside.plus)]
            index.inside.minus <- index.inside.minus[!is.na(index.inside.minus)]
            index <- c(index.inside.plus, index.inside.minus)
            index <- index[index>0 & index<=length(bins)] # index could cross chromosome boundaries, but we risk it
            if ('peaks' %in% what) binstates.inside[as.character(interval),] <- colMeans(binstates[index,])
            if ('combinations' %in% what) {
                fold <- table(combinations[index]) / tcombinations / length(annotation) * length(bins) # fold enrichment
#                 fold <- table(combinations[index]) / length(annotation) # percentage enrichment
                fold[is.na(fold)] <- 0
                combinations.inside[as.character(interval),] <- fold
            }
            if ('counts' %in% what) counts.inside[as.character(interval),] <- colMeans(counts[index,])
        }
        if ('peaks' %in% what) {
            enrich$peaks$inside <- binstates.inside
        }
        if ('combinations' %in% what) {
            enrich$combinations$inside <- combinations.inside
        }
        if ('counts' %in% what) {
            enrich$counts$inside <- counts.inside
        }
        stopTimedMessage(ptm)
    }

    ### 10000 bp before annotation ###
    if ('start' %in% region) {
        ptm <- startTimedMessage("Enrichment ",bp.around.annotation,"bp before annotations")
        # Get bins that overlap the start of annotation
        index.start.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="first")
        index.start.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="last")
        index.start.plus <- index.start.plus[!is.na(index.start.plus)]
        index.start.minus <- index.start.minus[!is.na(index.start.minus)]
        # Occurrences at every bin position relative to feature
        if ('peaks' %in% what) binstates.start <- array(dim=c(length(-lag:lag), length(info.dedup$ID)), dimnames=list(lag=-lag:lag, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.start <- array(dim=c(length(-lag:lag), length(levels(bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.start <- array(dim=c(length(-lag:lag), length(info$ID)), dimnames=list(lag=-lag:lag, track=info$ID))
        for (ilag in -lag:lag) {
            index <- c(index.start.plus+ilag, index.start.minus-ilag)
            index <- index[index>0 & index<=length(bins)]
            if ('peaks' %in% what) binstates.start[as.character(ilag),] <- colMeans(binstates[index,])
            if ('combinations' %in% what) {
                fold <- table(combinations[index]) / tcombinations / length(annotation) * length(bins) # fold enrichment
#                 fold <- table(combinations[index]) / length(annotation) # percentage enrichment
                fold[is.na(fold)] <- 0
                combinations.start[as.character(ilag),] <- fold
            }
            if ('counts' %in% what) counts.start[as.character(ilag),] <- colMeans(counts[index,])
        }
        if ('peaks' %in% what) {
            rownames(binstates.start) <- as.numeric(rownames(binstates.start)) * binsize
            enrich$peaks$start <- binstates.start
        }
        if ('combinations' %in% what) {
            rownames(combinations.start) <- as.numeric(rownames(combinations.start)) * binsize
            enrich$combinations$start <- combinations.start
        }
        if ('counts' %in% what) {
            rownames(counts.start) <- as.numeric(rownames(counts.start)) * binsize
            enrich$counts$start <- counts.start
        }
        stopTimedMessage(ptm)
    }

    ### 10000 bp after annotation ###
    if ('end' %in% region) {
        ptm <- startTimedMessage("Enrichment ",bp.around.annotation,"bp after annotations")
        # Get bins that overlap the end of annotation
        index.end.plus <- findOverlaps(annotation[strand(annotation)=='+' | strand(annotation)=='*'], bins, select="last")
        index.end.minus <- findOverlaps(annotation[strand(annotation)=='-'], bins, select="first")
        index.end.plus <- index.end.plus[!is.na(index.end.plus)]
        index.end.minus <- index.end.minus[!is.na(index.end.minus)]
        # Occurrences at every bin position relative to feature
        if ('peaks' %in% what) binstates.end <- array(dim=c(length(-lag:lag), length(info.dedup$ID)), dimnames=list(lag=-lag:lag, track=info.dedup$ID))
        if ('combinations' %in% what) combinations.end <- array(dim=c(length(-lag:lag), length(levels(bins$combination))), dimnames=list(lag=-lag:lag, combination=levels(bins$combination)))
        if ('counts' %in% what) counts.end <- array(dim=c(length(-lag:lag), length(info$ID)), dimnames=list(lag=-lag:lag, track=info$ID))
        for (ilag in -lag:lag) {
            index <- c(index.end.plus+ilag, index.end.minus-ilag)
            index <- index[index>0 & index<=length(bins)]
            if ('peaks' %in% what) binstates.end[as.character(ilag),] <- colMeans(binstates[index,])
            if ('combinations' %in% what) {
                fold <- table(combinations[index]) / tcombinations / length(annotation) * length(bins) # fold enrichment
#                 fold <- table(combinations[index]) / length(annotation) # percentage enrichment
                fold[is.na(fold)] <- 0
                combinations.end[as.character(ilag),] <- fold
            }
            if ('counts' %in% what) counts.end[as.character(ilag),] <- colMeans(counts[index,])
        }
        if ('peaks' %in% what) {
            rownames(binstates.end) <- as.numeric(rownames(binstates.end)) * binsize
            enrich$peaks$end <- binstates.end
        }
        if ('combinations' %in% what) {
            rownames(combinations.end) <- as.numeric(rownames(combinations.end)) * binsize
            enrich$combinations$end <- combinations.end
        }
        if ('counts' %in% what) {
            rownames(counts.end) <- as.numeric(rownames(counts.end)) * binsize
            enrich$counts$end <- counts.end
        }
        stopTimedMessage(ptm)
    }

    return(enrich[[what]])

}



