correlation.analysis <- function(multi.hmm, IDs=NULL, plot=TRUE) {

	## Intercept user input
	if (check.multivariate.model(multi.hmm)!=0) {
		message("Loading multivariate HMM from file ...", appendLF=F)
		multi.hmm <- get(load(multi.hmm))
		message(" done")
		if (check.multivariate.model(multi.hmm)!=0) stop("argument 'multi.hmm' expects a multivariate hmm object or a file that contains a multivariate hmm (type ?multi.hmm for help)")
	}
	if (is.null(IDs)) {
		IDs <- multi.hmm$IDs
	} else {
		if (length(IDs) != multi.hmm$num.modifications) {
			stop("argument 'IDs' needs one element for each mark in 'multi.hmm' (", multi.hmm$num.modifications, ")")
		}
	}

	## Calculate the correlation matrix on the combinatorial states
	binary.states <- dec2bin(as.integer(as.character(multi.hmm$states)))
	colnames(binary.states) <- IDs
	cor.matrix <- cor(binary.states)

	## Plot the correlation matrix
	if (!plot) {
		return(cor.matrix)
	} else {
		library(ggplot2)
		library(reshape2)
		mcor <- melt(cor.matrix, varnames=c("x","y"), value.name="correlation")
		ggplt <- ggplot() + geom_tile(data=mcor, mapping=aes_string(x='x', y='y', fill='correlation')) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient2()
		return(ggplt)
	}

}
