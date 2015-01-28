scanBinsizes = function(files.binned, outputfolder, chromosomes="chr10", eps=0.01, max.iter=100, max.time=300, repetitions=3, plot.progress=FALSE) {

	## Load libraries
	library(ggplot2)

	## Create outputfolder if not existent
	if (!file.exists(outputfolder)) dir.create(outputfolder)

	## Do univariates (real)
	path.uni = file.path(outputfolder, "results_univariate")
	if (!file.exists(path.uni)) dir.create(path.uni)

	for (binfile in files.binned) {
		unifile = file.path(path.uni, paste0("univariate_",basename(binfile)))
		if (!file.exists(unifile)) {
			binned.data <- get(load(binfile))
			if (!is.null(chromosomes)) {
				binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
			}
			unimodel = callPeaksUnivariate(binned.data, ID=basename(binfile), eps=eps, max.iter=max.iter, max.time=max.time)
			save(unimodel, file=unifile)
		}
	}

	performance = NULL
	for (irep in 1:repetitions) {

		## Simulate data from univariates
		cat("#--- Simulate data ---#\n")
		path.sim.bin = file.path(outputfolder, paste0("simulated_rep_",irep,"_binned"))
		if (!file.exists(path.sim.bin)) dir.create(path.sim.bin)
		simfile = file.path(path.sim.bin, paste0("simulated.RData"))
		if (!file.exists(simfile)) {
			files2sim = list.files(path.uni, full.names=TRUE)
			sim.data.list = list()
			for (binfile in files.binned) {
				cat(binfile,"\n")
				unifile = file.path(path.uni, paste0("univariate_",basename(binfile)))
				unimodel = get(load(unifile))
				sim.data = simulateUnivariate(unimodel$bins, unimodel$transitionProbs, unimodel$distributions)
				sim.data.list[[basename(binfile)]] = sim.data
				# Save to separate binfiles
				binsize = width(sim.data$bins)[1]
				simulated.binned.data <- sim.data$bins
				save(simulated.binned.data, file=file.path(path.sim.bin, paste0("simulated_",basename(binfile))))
			}
			# Save with states
			save(sim.data.list, file=simfile)
		} else {
			sim.data.list = get(load(simfile))
		}
		# Get simulated original states
		cat("Getting simulated original states...")
		sim.states.list = lapply(lapply(sim.data.list, "[[", 'bins'), function(gr) { return(gr$state) } )
		for (i1 in 1:length(sim.states.list)) {
			sim.states.list[[i1]] = state.labels[sim.states.list[[i1]]]
		}
		cat(" done\n")

		## Univariate (simulated)
		path.sim.uni = file.path(outputfolder, paste0("simulated_rep_",irep,"_results_univariate"))
		if (!file.exists(path.sim.uni)) dir.create(path.sim.uni)
		sim.binfiles = list.files(path.sim.bin, full.names=TRUE, pattern="simulated_")

		# Go through binsizes
		uni.states.list = list()
		binsizes <- vector()
		for (binfile in sim.binfiles) {
			unifile = file.path(path.sim.uni, paste0("univariate_",basename(binfile)))
			if (!file.exists(unifile)) {
				binned.data <- get(load(binfile))
				if (!is.null(chromosomes)) {
					binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
				}
				unimodel = callPeaksUnivariate(binned.data, ID=basename(binfile), eps=eps, max.iter=max.iter, max.time=max.time)
				save(unimodel, file=unifile)
			} else {
				unimodel <- get(load(unifile))
			}
			# Get univariate predicted states
			cat("Getting predicted univariate states...")
			uni.states.list[[i1]] = unimodel$bins$state
			binsizes[i1] <- width(unimodel$bins)[1]
			cat(" done\n")

			# Performance
			cat("Calculating performance...")
			mask = sim.states.list[[i1]] != uni.states.list[[i1]]
			miscalls = length(which(mask)) / length(sim.states.list[[i1]])
			performance = rbind(performance, data.frame(binsize=binsizes[i1], miscalls, state1weight=unimodel$weights[3]))
			cat(" done\n")

			## Plot miscalls
			if (plot.progress) {
				ggplt = ggplot() + theme_bw() + geom_boxplot(data=performance, aes(x=as.factor(binsize), y=miscalls)) + scale_x_discrete(limits=rev(levels(performance$binsize))) + labs(title=basename(outputfolder)) + xlab('binsize') + ylab('fraction of miscalls')
				print(ggplt)
			}
		}
	}

	rownames(performance) <- NULL
	ggplt = ggplot() + theme_bw() + geom_boxplot(data=performance, aes(x=as.factor(binsize), y=miscalls)) + scale_x_discrete(limits=rev(levels(performance$binsize))) + labs(title=basename(outputfolder)) + xlab('binsize') + ylab('fraction of miscalls')
	return(list(performance, ggplt))

}
