find.optimal.binsize = function(bedfile, reference.genome.file, outputfolder, binsizes=c(seq(10000,4000,-2000),seq(3600,400,-400),200,100), chromosome="chr10", repetitions=3, plot.progress=FALSE) {

	## Create outputfolder if not existent
	if (!file.exists(outputfolder)) dir.create(outputfolder)

	## Bin the data
	path.bin = file.path(outputfolder,"real_binned_data")
	if (length(list.files(path.bin)) < length(binsizes)) {
		bin.BED(bedfile,reference.genome.file,outputfolder=path.bin,binsizes=binsizes,chromosome=chromosome)
	}

	## Univariate (real)
	path.uni = file.path(outputfolder,"real_uniresults")
	if (!file.exists(path.uni)) dir.create(path.uni)
	files2uni = rev(mixedsort(list.files(path.bin, full=TRUE)))

	for (i1 in 1:length(files2uni)) {
		binfile = files2uni[i1]
		unifile = paste(path.uni,.Platform$file.sep,"uniresult_",basename(binfile), sep="")
		if (!file.exists(unifile)) {
			binned.data = get(load(binfile))
			unimodel = univariate.from.binned.data(binned.data, eps=0.01)
			unimodel$binfile = basename(binfile)
			save(unimodel, file=unifile)
		}
	}

	performance = NULL
	for (irep in 1:repetitions) {

		## Simulate data from univariates
		cat("#--- Simulate data ---#\n")
		path.sim = paste(outputfolder,.Platform$file.sep,"sim_binned_data_",irep, sep="")
		if (!file.exists(path.sim)) dir.create(path.sim)
		simfile = paste(path.sim,.Platform$file.sep,"sim_data.rda", sep="")
		if (!file.exists(simfile)) {
			files2sim = rev(mixedsort(list.files(path.uni, full=TRUE)))
			sim.data.list = list()
			for (i1 in 1:length(files2sim)) {
				cat("binsize =",binsizes[i1],"\n")
				unifile = files2sim[i1]
				unimodel = get(load(unifile))
				sim.data = simulate.univariate(unimodel$coordinates, unimodel$A, unimodel$distributions)
				sim.data.list[[i1]] = sim.data
				# Save to separate binfiles
				binsize = sim.data$coordinates$start[2]
				binned.data = data.frame(sim.data$coordinates, sim.data$reads)
				save(binned.data, file=paste(path.sim,.Platform$file.sep,"simulated_binsize_",binsize,"_",chromosome,".RData", sep=""))
			}
			# Save with states
			save(sim.data.list, file=simfile)
		} else {
			sim.data.list = get(load(simfile))
		}
		# Get simulated original states
		cat("Getting simulated original states...")
		sim.states.list = lapply(sim.data.list, "[[", 'states')
		for (i1 in 1:length(sim.states.list)) {
			sim.states.list[[i1]] = ifelse(sim.states.list[[i1]]==3, 1, 0)
		}
		cat(" done\n")

		## Univariate (simulated)
		path.bin = paste(outputfolder,.Platform$file.sep,"sim_binned_data_",irep, sep="")
		path.uni = paste(outputfolder,.Platform$file.sep,"sim_uniresults_",irep, sep="")
		if (!file.exists(path.uni)) dir.create(path.uni)
		files2uni = rev(mixedsort(list.files(path.bin, full=TRUE, pattern="simulated")))

		# Go through binsizes
		uni.states.list = list()
		for (i1 in 1:length(files2uni)) {
			binfile = files2uni[i1]
			unifile = paste(path.uni,.Platform$file.sep,"uniresult_",basename(binfile), sep="")
			if (!file.exists(unifile)) {
				binned.data = get(load(binfile))
				unimodel = univariate.from.binned.data(binned.data, eps=0.1)
				save(unimodel, file=unifile)
			} else {
				unimodel = get(load(unifile))
			}
			# Get univariate predicted states
			cat("Getting predicted univariate states...")
			uni.states.list[[i1]] = get.states(unimodel)
			cat(" done\n")

			# Performance
			cat("Calculating performance...")
			mask = sim.states.list[[i1]] != uni.states.list[[i1]]
			miscalls = length(which(mask)) / length(sim.states.list[[i1]])
			performance = rbind(performance, data.frame(binsize=as.factor(binsizes[i1]), miscalls, state1weight=unimodel$softweights[3]))
			cat(" done\n")

			## Plot miscalls
			if (plot.progress) {
				ggplt = ggplot() + theme_bw() + geom_boxplot(data=performance, aes(x=binsize, y=miscalls)) + scale_x_discrete(limits=rev(levels(performance$binsize))) + labs(title=basename(outputfolder))
# + scale_x_discrete(limits=as.character(seq(min(binsizes),max(binsizes),min(abs(diff(binsizes))))))
				print(ggplt)
			}
		}
	}

	return(performance)

}
