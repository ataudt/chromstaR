gwasCatalog2GRanges <- function(gwas.catalog.file) {
	
	data.unf <- read.csv(gwas.catalog.file, header=TRUE, sep="\t")
	# Filter NA rows
	data <- data.unf[!is.na(data.unf[,"Chr_id"]) & !is.na(data.unf[,"Chr_pos"]) , ]
	gr <- GenomicRanges::GRanges(seqnames = paste0("chr",data[,"Chr_id"]),
																ranges = IRanges(start=data[,"Chr_pos"], end=data[,"Chr_pos"]),
																strand = Rle(strand("*"), nrow(data)),
																reported.gene = data[,"Reported.Gene.s."],
																intergenic = data[,"Intergenic"]
																)

	return(gr)


}

enrichment.from.expression <- function(multi.hmm, bamfile, bamindex=bamfile, per.mark=FALSE) {

	## Intercept user input
	if (check.multivariate.model(multi.hmm)!=0) {
		cat("Loading multivariate HMM from file ...")
		multi.hmm <- get(load(multi.hmm))
		cat(" done\n")
		if (check.multivariate.model(multi.hmm)!=0) stop("argument 'multi.hmm' expects a multivariate hmm object or a file that contains a multivariate hmm (type ?multi.hmm for help)")
	}

	## Convert hmm to GRanges
	gr <- hmm2GRanges(multi.hmm, reduce=FALSE)

	## Read in BAM file
	align <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(what=c("pos"),which=range(gr)))

	## Create annotation object with same coordinates as HMM
	anno <- gr
	mcols(anno) <- NULL

	### Categorize the read count
	# To get a good categorization, split the bins into unexpressed (zero-inflation + unmodified) and expressed (modified) bins with our HMM
	mcols(anno)$reads <- countOverlaps(anno, align)
	expr.hmm <- univariate.from.binned.data(anno, ID="expr", max.iter=100, eps=1, output.if.not.converged=TRUE)
	## Set the categories based on the expressed bins
	expr.breaks <- quantile(expr.hmm$reads[expr.hmm$states=='modified'], 0:10/11)
	# Insert interval [0,1) into breaks if not there
	expr.breaks <- unique(sort(c(0,1,expr.breaks)))
	expr.category <- findInterval(mcols(anno)$reads, expr.breaks)

	labels <- paste0("expr.", expr.breaks[expr.category], "-", expr.breaks[expr.category+1])
	mcols(anno)$feature <- factor(labels, levels=mixedsort(unique(labels)))

	## Do the enrichment analysis
	cat("Doing enrichment analysis ...")
	if (per.mark) {
		enrichment.table <- NULL
		gr.mark <- gr
		mcols(gr.mark) <- NULL
		# Convert combinatorial states (factors) to binary representation
		binary.states <- dec2bin(as.integer(as.character(mcols(gr)$state)))
		colnames(binary.states) <- multi.hmm$IDs.univariate
		# Do enrichment analysis for every mark
		for (col in 1:ncol(binary.states)) {
			mcols(gr.mark)$state <- as.integer(binary.states[,col])
			enrichment.table.mark <- enrichment.analysis(gr.mark, anno)
			# Select only 'modified' row
			enrichment.table.mark <- enrichment.table.mark[2,]
			# Set name properly
			enrichment.table.mark$state <- colnames(binary.states)[col]
			# Combine to result
			enrichment.table <- rbind(enrichment.table, enrichment.table.mark)
		}
	} else {
		enrichment.table <- enrichment.analysis(gr, anno)
	}
	cat(" done\n")

	## Return results
	return(enrichment.table)


}

enrichment.from.annotation <- function(multi.hmm, annotation.file.gtf, per.mark=FALSE) {

	## Intercept user input
	if (check.multivariate.model(multi.hmm)!=0) {
		cat("Loading multivariate HMM from file ...")
		multi.hmm <- get(load(multi.hmm))
		cat(" done\n")
		if (check.multivariate.model(multi.hmm)!=0) stop("argument 'multi.hmm' expects a multivariate hmm object or a file that contains a multivariate hmm (type ?multi.hmm for help)")
	}

	## Load annotation file
	library(rtracklayer)
	cat("Loading annotation file ...")
	anno <- rtracklayer::import(annotation.file.gtf, format="gtf")
	cat(" done\n")

	## Convert hmm to GRanges
	gr <- hmm2GRanges(multi.hmm, reduce=FALSE)

	## Do the enrichment analysis
	cat("Doing enrichment analysis ...")
	mcols(anno)$feature <- mcols(anno)$type
	if (per.mark) {
		enrichment.table <- NULL
		gr.mark <- gr
		mcols(gr.mark) <- NULL
		# Convert combinatorial states (factors) to binary representation
		binary.states <- dec2bin(as.integer(as.character(mcols(gr)$state)))
		colnames(binary.states) <- multi.hmm$IDs.univariate
		# Do enrichment analysis for every mark
		for (col in 1:ncol(binary.states)) {
			mcols(gr.mark)$state <- as.integer(binary.states[,col])
			enrichment.table.mark <- enrichment.analysis(gr.mark, anno)
			# Select only 'modified' row
			enrichment.table.mark <- enrichment.table.mark[2,]
			# Set name properly
			enrichment.table.mark$state <- colnames(binary.states)[col]
			# Combine to result
			enrichment.table <- rbind(enrichment.table, enrichment.table.mark)
		}
	} else {
		enrichment.table <- enrichment.analysis(gr, anno)
	}
	cat(" done\n")

	## Return results
	return(enrichment.table)

}

enrichment.analysis <- function(granges.states.per.bin, granges.annotation) {

	# Rename variables for easier access
	gr.states <- granges.states.per.bin
	gr.anno <- granges.annotation

	# Split gr.states into list by states
	grl.states <- split(gr.states, mcols(gr.states)$state)

	# Total number of bins
	num.bins <- length(gr.states)

	# Features in gr.annotation
	features <- levels(mcols(gr.anno)$feature)

	# P-values for overlaps into each combinatorial state
	overlaps.per.state <- NULL
	num.feature <- NULL
	num.nonfeature <- NULL
	num.bins.per.state <- unlist(lapply(grl.states, length))
	p.depleted <- NULL
	p.enriched <- NULL
	fold.enrichment <- NULL
	for (feature in features) {
		i1 <- which(feature==features)
		mask <- mcols(gr.anno)$feature == feature
		overlaps.per.state[[i1]] <- countOverlaps(grl.states, gr.anno[mask])
		num.feature[[i1]] <- sum(overlaps.per.state[[i1]])
		num.nonfeature[[i1]] <- num.bins - num.feature[[i1]]
		# p-values
		p.depleted[[i1]] <- phyper(overlaps.per.state[[i1]], num.feature[[i1]], num.nonfeature[[i1]], num.bins.per.state)
		p.enriched[[i1]] <- 1 - phyper(overlaps.per.state[[i1]]-1, num.feature[[i1]], num.nonfeature[[i1]], num.bins.per.state)
		# fold enrichment
		fold.enrichment[[i1]] <- log( base=2, overlaps.per.state[[i1]] / num.feature[[i1]] / num.bins.per.state * num.bins )
	}
	# Multiple testing correction
	p.depleted <- matrix(p.adjust(unlist(p.depleted), method="bonferroni"), ncol=length(features))
	colnames(p.depleted) <- features
	rownames(p.depleted) <- names(grl.states)
	p.enriched <- matrix(p.adjust(unlist(p.enriched), method="bonferroni"), ncol=length(features))
	colnames(p.enriched) <- features
	rownames(p.enriched) <- names(grl.states)

	## Make return data.frame
	enrich <- data.frame(
		state = names(grl.states),
		num.bins = num.bins.per.state,
		perc.bins = num.bins.per.state / num.bins * 100,
		overlaps.per.state,
		lapply(lapply(overlaps.per.state, "/", num.bins.per.state), "*", 100),
		fold.enrichment,
		p.depleted,
		p.enriched
	)
	names(enrich) <- c(
		"state",
		"num.bins",
		"%.bins",
		paste("num", features, sep="."),
		paste("%", features, sep="."),
		paste(features, "fold_enrich", sep="."),
		paste(features, "p_dep", sep="."),
		paste(features, "p_enr", sep=".")
	)
	# Reorder columns
# 	lf = length(features)
# 	reorder = c(1:3, as.vector( rbind( (4):(4+1*lf-1),(4+1*lf):(4+2*lf-1),(4+2*lf):(4+3*lf-1),(4+3*lf):(4+4*lf-1) ) ) )
# 	enrich = enrich[reorder]
	enrich <- as.data.frame(enrich)
	enrich <- enrich[order(enrich$num.bins, decreasing=TRUE),]
	enrich <- apply(enrich, 2, as.double)
	enrich <- as.data.frame(enrich)

	return(enrich)

}




# enrichment.analysis = function(annotation.table, comb.states) {
# 
# 	## Check user input
# 	if (nrow(annotation.table)!=length(comb.states)) {
# 		stop("both arguments must have the same length")
# 	}
# 
# 	## Insert artificial column to get 'bincount' later
# 	numbins = length(comb.states) # total number of bins
# 	# Delete coordinate columns and add bincount
# 	drops = c("chrom","start","end")
# 	bincount = rep(1,nrow(annotation.table)) # contains the number of bins in each states later
# 	annotation.table = data.frame(annotation.table[ , !(names(annotation.table) %in% drops)], bincount = bincount)
# 
# 	## Convert numbers to logicals (a feature is either present or not at a given position)
# 	cat("convert to logical...      \r")
# 	annotations2aggregate = as.data.frame(lapply(annotation.table, as.logical))
# 	## Aggregate
# 	cat("aggregating...             \r")
# 	aggregated = aggregate(annotations2aggregate, by=list(comb.state = comb.states), sum)
# 	ind_aggregated = 2:(ncol(aggregated)-1) # first columns now holds the comb.states, last is bincount
# 
# 	## Calculate p-values
# 	cat("calculate p-values...      \r")
# 	p.depleted = NULL
# 	p.enriched = NULL
# 	for (icol in ind_aggregated) {
# 		x = aggregated[,icol] # number of overlaps
# 		m = sum(x) # number of feature in column
# 		n = numbins - m # number of not-feature in column
# 		k = aggregated[,"bincount"] # number of combinatorial states
# 		p.depleted[[length(p.depleted)+1]] = phyper(x,m,n,k)
# 		p.enriched[[length(p.enriched)+1]] = 1 - phyper(x-1,m,n,k)
# 	}
# 	p.depleted = matrix(p.adjust(unlist(p.depleted), method="holm"), ncol=length(ind_aggregated))
# 	p.enriched = matrix(p.adjust(unlist(p.enriched), method="holm"), ncol=length(ind_aggregated))
# 
# 	## Make return data.frame
# 	cat("concatenate...             \r")
# 	out = data.frame(
# 		comb.state = aggregated[,"comb.state"],
# 		num.bins = aggregated[,"bincount"],
# 		perc.bins = aggregated[,"bincount"] / numbins * 100,
# 		aggregated[,ind_aggregated],
# 		aggregated[,ind_aggregated] / aggregated[,"bincount"] * 100,
# 		p.depleted,
# 		p.enriched
# 	)
# 	names(out) = c(
# 		"comb.state",
# 		"bins.num",
# 		"bins.%.data",
# 		paste(names(aggregated[,ind_aggregated]), "num", sep="."),
# 		paste(names(aggregated[,ind_aggregated]), "%", "state", sep="."),
# 		paste(names(aggregated[,ind_aggregated]), "p_depleted", sep="."),
# 		paste(names(aggregated[,ind_aggregated]), "p_enriched", sep=".")
# 	)
# 	# Reorder columns
# 	lia = length(ind_aggregated)
# 	reorder = c(1:3, as.vector( rbind( (4):(4+1*lia-1),(4+1*lia):(4+2*lia-1),(4+2*lia):(4+3*lia-1),(4+3*lia):(4+4*lia-1) ) ) )
# 	out = out[reorder]
# 	cat("                           \r")
# 
# 	return(out)
# 
# }
