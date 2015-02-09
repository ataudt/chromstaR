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
			message()
			message.underlined(paste0('Univariate for ',basename(binfile),':'), line='=')
			binned.data <- get(load(binfile))
			if (!is.null(chromosomes)) {
				binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
			}
			unimodel = callPeaksUnivariate(binned.data, ID=basename(binfile), eps=eps, max.iter=max.iter, max.time=max.time)
			save(unimodel, file=unifile)
		}
	}

	message()
	performance = NULL
	for (irep in 1:repetitions) {
		message.underlined(paste0('Simulation ',irep,':'), line='=')

		## Simulate data from univariates
		path.sim.bin = file.path(outputfolder, paste0("simulated_rep_",irep,"_binned"))
		if (!file.exists(path.sim.bin)) dir.create(path.sim.bin)
		simfile = file.path(path.sim.bin, paste0("simulated.RData"))
		if (!file.exists(simfile)) {
			message('  simulating data ...', appendLF=F); ptm <- proc.time()
			files2sim = list.files(path.uni, full.names=TRUE)
			sim.data.list = list()
			for (binfile in files.binned) {
				unifile = file.path(path.uni, paste0("univariate_",basename(binfile)))
				unimodel = get(load(unifile))
				sim.data = suppressMessages( simulateUnivariate(unimodel$bins, unimodel$transitionProbs, unimodel$distributions) )
				sim.data.list[[basename(binfile)]] = sim.data
				# Save to separate binfiles
				binsize = width(sim.data$bins)[1]
				simulated.binned.data <- sim.data$bins
				sim.binfile <- file.path(path.sim.bin, paste0("simulated_",basename(binfile)))
				save(simulated.binned.data, file=sim.binfile)
			}
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
			# Save with states
			message('  saving to file ...', appendLF=F); ptm <- proc.time()
			save(sim.data.list, file=simfile)
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		} else {
			message('  loading simulated data ...', appendLF=F); ptm <- proc.time()
			sim.data.list = get(load(simfile))
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		}
		# Get simulated original states
		message("  getting simulated original states ...", appendLF=F); ptm <- proc.time()
		sim.states.list = lapply(lapply(sim.data.list, "[[", 'bins'), function(gr) { return(gr$state) } )
		for (binfile in files.binned) {
			sim.states.list[[basename(binfile)]] = state.labels[sim.states.list[[basename(binfile)]]]
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		## Univariate (simulated)
		path.sim.uni = file.path(outputfolder, paste0("simulated_rep_",irep,"_results_univariate"))
		if (!file.exists(path.sim.uni)) dir.create(path.sim.uni)

		# Go through binsizes (=binfiles)
		uni.states.list = list()
		binsizes <- vector()
		for (binfile in files.binned) {
			sim.binfile <- file.path(path.sim.bin, paste0('simulated_',basename(binfile)))
			sim.unifile = file.path(path.sim.uni, paste0("univariate_simulated_",basename(binfile)))
			if (!file.exists(sim.unifile)) {
				message()
				message.underlined(paste0('Univariate for simulated ',basename(binfile),':'), line='-')
				binned.data <- get(load(sim.binfile))
				if (!is.null(chromosomes)) {
					binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
				}
				unimodel = callPeaksUnivariate(binned.data, ID=basename(sim.binfile), eps=eps, max.iter=max.iter, max.time=max.time)
				save(unimodel, file=sim.unifile)
			} else {
				unimodel <- get(load(sim.unifile))
			}
			# Get univariate predicted states
			message("  getting predicted univariate states...", appendLF=F); ptm <- proc.time()
			uni.states.list[[basename(binfile)]] = unimodel$bins$state
			binsizes[basename(binfile)] <- width(unimodel$bins)[1]
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

			# Performance
			message("  calculating performance...", appendLF=F); ptm <- proc.time()
			mask = sim.states.list[[basename(binfile)]] != uni.states.list[[basename(binfile)]]
			miscalls = length(which(mask)) / length(sim.states.list[[basename(binfile)]])
			performance = rbind(performance, data.frame(simulation=irep, binsize=binsizes[basename(binfile)], miscalls))
			time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

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
