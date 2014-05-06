enrichment.analysis = function(annotation.table, comb.states) {

	## Check user input
	if (nrow(annotation.table)!=length(comb.states)) {
		stop("both arguments must have the same length")
	}

	## Insert artificial column to get 'bincount' later
	numbins = length(comb.states) # total number of bins
	# Delete coordinate columns and add bincount
	drops = c("chrom","start","end")
	bincount = rep(1,nrow(annotation.table)) # contains the number of bins in each states later
	annotation.table = data.frame(annotation.table[ , !(names(annotation.table) %in% drops)], bincount = bincount)

	## Convert numbers to logicals (a feature is either present or not at a given position)
	cat("convert to logical...      \r")
	annotations2aggregate = as.data.frame(lapply(annotation.table, as.logical))
	## Aggregate
	cat("aggregating...             \r")
	aggregated = aggregate(annotations2aggregate, by=list(comb.state = comb.states), sum)
	ind_aggregated = 2:(ncol(aggregated)-1) # first columns now holds the comb.states, last is bincount

	## Calculate p-values
	cat("calculate p-values...      \r")
	p.depleted = NULL
	p.enriched = NULL
	for (icol in ind_aggregated) {
		x = aggregated[,icol] # number of overlaps
		m = sum(x) # number of feature in column
		n = numbins - m # number of not-feature in column
		k = aggregated[,"bincount"] # number of combinatorial states
		p.depleted[[length(p.depleted)+1]] = phyper(x,m,n,k)
		p.enriched[[length(p.enriched)+1]] = 1 - phyper(x-1,m,n,k)
	}
	p.depleted = matrix(p.adjust(unlist(p.depleted), method="holm"), ncol=length(ind_aggregated))
	p.enriched = matrix(p.adjust(unlist(p.enriched), method="holm"), ncol=length(ind_aggregated))

	## Make return data.frame
	cat("concatenate...             \r")
	out = data.frame(
		comb.state = aggregated[,"comb.state"],
		num.bins = aggregated[,"bincount"],
		perc.bins = aggregated[,"bincount"] / numbins * 100,
		aggregated[,ind_aggregated],
		aggregated[,ind_aggregated] / aggregated[,"bincount"] * 100,
		p.depleted,
		p.enriched
	)
	names(out) = c(
		"comb.state",
		"bins.num",
		"bins.%.data",
		paste(names(aggregated[,ind_aggregated]), "num", sep="."),
		paste(names(aggregated[,ind_aggregated]), "%", "state", sep="."),
		paste(names(aggregated[,ind_aggregated]), "p_depleted", sep="."),
		paste(names(aggregated[,ind_aggregated]), "p_enriched", sep=".")
	)
	# Reorder columns
	lia = length(ind_aggregated)
	reorder = c(1:3, as.vector( rbind( (4):(4+1*lia-1),(4+1*lia):(4+2*lia-1),(4+2*lia):(4+3*lia-1),(4+3*lia):(4+4*lia-1) ) ) )
	out = out[reorder]
	cat("                           \r")

	return(out)

}
