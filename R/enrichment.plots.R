# plotFractionOfSegmentation <- function(multi.hmm) {
# 
#     # Convert to GRanges if not already GRanges
#     if (class(multi.hmm)=="GRanges") {
#         gr <- multi.hmm
#     } else if (class(multi.hmm)==class.multivariate.hmm) {
#         gr <- multi.hmm$segments
#     } else {
#         stop("argument 'multi.hmm' expects either a multivariate hmm object (type ?multi.hmm for help) or a derived GRanges object")
#     }
#     grl <- GenomicRanges::split(gr, GenomicRanges::mcols(gr)$state)
#     # Extract the information
#     num.segments.per.state <- lapply(grl, length)
#     bp.per.state <- lapply(grl, function(x) { sum(as.numeric(width(ranges(x)))) } )
#     df <- data.frame(segments = unlist(num.segments.per.state), 
#                                         bases = unlist(bp.per.state))
#     df <- sweep(df, 2, colSums(df), "/")
#     df$state <- rownames(df)
# 
#     # Plot
#     dfplot <- melt(df, id.vars='state', variable.name='category', value.name='fraction')
# #     dfplot$state <- factor(dfplot$state, levels=df$state[rev(order(df$bases))])
#     dfplot$state <- factor(dfplot$state, levels=levels(mcols(gr)$state))
#     ggplt <- ggplot(data=dfplot) + theme_bw() + geom_bar(aes_string(x='state', y='fraction', fill='category'), position='dodge', stat='identity') + scale_fill_discrete(guide = guide_legend(reverse=TRUE, title=NULL)) + theme(legend.position='top') + ylab('fraction of segmentation') + xlab('combinatorial state') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#     
#     return(ggplt)
# 
# }
# 
# 
# plotEnrichCountHeatmapFoldEnrichment <- function(enrichment.table) {
# 
#     et <- enrichment.table
# 
#     ## Get columns with fold enrichment values
#     fe.names <- unlist(lapply(as.list(names(et)), function(x) { grep(x=x, pattern='fold_enrich', value=TRUE) } ))
#     fe <- et[,c('state',fe.names)]
#     names(fe) <- gsub(pattern='\\.fold_enrich', replacement='', x=names(fe))
#     df.feat <- melt(fe, id.vars='state', variable.name='feature', value.name='fold.enrichment')
# 
#     ## Get columns with p-values
#     p.dep.names <- unlist(lapply(as.list(names(et)), function(x) { grep(x=x, pattern='p_dep', value=TRUE) } ))
#     pv <- et[,c('state',p.dep.names)]
#     df.p.dep <- melt(pv, id.vars='state', variable.name='feature', value.name='p.value')
#     p.enr.names <- unlist(lapply(as.list(names(et)), function(x) { grep(x=x, pattern='p_enr', value=TRUE) } ))
#     pv <- et[,c('state',p.enr.names)]
#     df.p.enr <- melt(pv, id.vars='state', variable.name='feature', value.name='p.value')
#     # Combine
#     df.feat$significant <- rep('', nrow(df.feat))
#     df.feat$significant[df.p.dep$p.value < 0.05] <- 'd'
#     df.feat$significant[df.p.enr$p.value < 0.05] <- 'E'
# 
#     ## Heatmap
#     df.feat$state <- factor(df.feat$state, levels=fe$state)
#     df.feat$feature <- factor(df.feat$feature, levels=mixedsort(unique(df.feat$feature)))
#     ggplt <- ggplot(data=df.feat) + geom_tile(aes_string(y='feature', x='state', fill='fold.enrichment')) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient2(low='blue', mid='yellow', high='red', na.value='white') + geom_text(aes_string(y='feature', x='state', label='significant'), size=3)
# 
#     return(ggplt)
# }
