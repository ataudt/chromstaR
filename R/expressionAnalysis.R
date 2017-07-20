#' Overlap with expression data
#'
#' Get the expression values that overlap with each combinatorial state.
#'
#' @param hmm A \code{\link{multiHMM}} or \code{\link{combinedMultiHMM}} object or file that contains such an object.
#' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range.
#' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
#' @param return.marks Set to \code{TRUE} if expression values for marks instead of combinations should be returned.
#' @return A \code{\link{ggplot2}} object if a \code{\link{multiHMM}} was given or a named list with \code{\link{ggplot2}} objects if a \code{\link{combinedMultiHMM}} was given.
#' @author Aaron Taudt
#' @seealso \code{\link{plotting}}
#' @importFrom IRanges subsetByOverlaps
#' @importFrom reshape2 melt
#' @export
#' @examples
#' ## Load an example multiHMM
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
#'
#' ## Obtain expression data
#'data(expression_lv)
#'head(expression_lv)
#'
#'## We need to get coordinates for each of the genes
#'library(biomaRt)
#'ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host='may2012.archive.ensembl.org',
#'                   dataset='rnorvegicus_gene_ensembl')
#'genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position',
#'                            'end_position', 'strand', 'external_gene_id',
#'                            'gene_biotype'),
#'               mart=ensembl)
#'expr <- merge(genes, expression_lv, by='ensembl_gene_id')
#'# Transform to GRanges
#'expression.SHR <- GRanges(seqnames=paste0('chr',expr$chromosome_name),
#'                          ranges=IRanges(start=expr$start, end=expr$end),
#'                          strand=expr$strand, name=expr$external_gene_id,
#'                          biotype=expr$gene_biotype,
#'                          expression=expr$expression_SHR)
#'# We apply an asinh transformation to reduce the effect of outliers
#'expression.SHR$expression <- asinh(expression.SHR$expression)
#'
#'## Plot
#'plotExpression(model, expression.SHR) +
#'  theme(axis.text.x=element_text(angle=0, hjust=0.5)) +
#'  ggtitle('Expression of genes overlapping combinatorial states')
#'plotExpression(model, expression.SHR, return.marks=TRUE) +
#'  ggtitle('Expression of marks overlapping combinatorial states')
#'
plotExpression <- function(hmm, expression, combinations=NULL, return.marks=FALSE) {
    
    hmm <- loadHmmsFromFiles(hmm, check.class=c(class.multivariate.hmm, class.combined.multivariate.hmm))[[1]]
    ## Variables
    bins <- hmm$bins
    if (class(hmm) == class.combined.multivariate.hmm) {
    } else if (class(hmm) == class.multivariate.hmm) {
        # Rename 'combination' to 'combination.' for coherence with combinedMultiHMM
        names(mcols(bins))[grep('combination', names(mcols(bins)))] <- paste0('combination.', unique(hmm$info$condition))
    }
    conditions <- sub('combination.', '', grep('combination', names(mcols(bins)), value=TRUE))
    if (is.null(combinations)) {
        comb.levels <- levels(mcols(bins)[,paste0('combination.', conditions[1])])
    } else {
        comb.levels <- combinations
    }
    marks <- unique(unlist(strsplit(gsub('\\]','',gsub('\\[','',comb.levels)),'\\+')))
    
    ggplts <- list()
    for (condition in conditions) {
        bins$combination <- mcols(bins)[,paste0('combination.', condition)]
        exprlist <- list()
        if (return.marks) {
            for (mark in marks) {
                mask <- grepl(paste0('\\<',mark,'\\>'),bins$combination)
                expr.mark <- IRanges::subsetByOverlaps(expression, bins[mask])
                exprlist[[mark]] <- expr.mark$expression
            }
        } else {
            for (comb.level in comb.levels) {
                mask <- bins$combination == comb.level
                expr.combstate <- IRanges::subsetByOverlaps(expression, bins[mask])
                exprlist[[comb.level]] <- expr.combstate$expression
            }
        }
        
        df <- reshape2::melt(exprlist)
        if (return.marks) {
            names(df) <- c('expression', 'mark')
            ggplt <- ggplot(df) + geom_boxplot(aes_string(x='mark', y='expression'))
        } else {
            names(df) <- c('expression', 'combination')
            ggplt <- ggplot(df) + geom_boxplot(aes_string(x='combination', y='expression'))
        }
        ggplt <- ggplt + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        ggplts[[condition]] <- ggplt
    }
    
    if (class(hmm) == class.multivariate.hmm) {
        return(ggplts[[1]])
    } else if (class(hmm) == class.combined.multivariate.hmm) {
        return(ggplts)
    }

}


# #' Mean expression at percentage overlap
# #'
# #' Get the average expression for each percentage of overlap of combinatorial state with feature.
# #'
# #' @param multi.hmm A \code{\link{multiHMM}} or a file that contains such an object.
# #' @param expression A \code{\link{GRanges}} object with metadata column 'expression', containing the expression value for each range of the feature.
# #' @param combinations A vector with combinations for which the expression overlap will be calculated. If \code{NULL} all combinations will be considered.
# #' @return A list with vectors of mean expression values per percentile for each combinatorial state. 
# #' @author Aaron Taudt
# #' @importFrom S4Vectors queryHits
# #' @export
# expressionAtPercentageOverlap <- function(multi.hmm, expression, combinations=NULL) {
# 
#     multi.hmm <- loadHmmsFromFiles(multi.hmm, check.class=class.multivariate.hmm)[[1]]
#     ## Variables
#     bins <- multi.hmm$bins
#     if (is.null(combinations)) {
#         comb.levels <- levels(bins$combination)
#     } else {
#         comb.levels <- combinations
#     }
#     nintervals <- 100
#     
#     expression.means <- array(NA, dim=c(nintervals+1, length(comb.levels), 2), dimnames=list(percentage=0:nintervals, combination=comb.levels, value=c('expression','weight')))
#     for (icomb in 1:length(comb.levels)) {
#         mask <- bins$combination == comb.levels[icomb]
#         ind <- findOverlaps(expression, bins[mask])
#         rle <- rle(S4Vectors::queryHits(ind))
#         expression$num.bins <- 0
#         expression$num.bins[rle$values] <- rle$lengths
#         expression$genewidth <- width(expression)
#         expression$percentage <- round(expression$num.bins*1000 / expression$genewidth * nintervals) # Normalize to genewidth
#         expression$percentage[expression$percentage>=nintervals] <- nintervals
#         splt <- split(expression$expression, expression$percentage)
#         tab <- sapply(splt, mean, na.rm=TRUE)
#         expression.means[names(tab),icomb,'expression'] <- tab #select by icomb instead of name because of potential '' states
#         expression.means[names(tab),icomb,'weight'] <- sapply(splt, length)
#     }
#     return(expression.means)
# 
# }
