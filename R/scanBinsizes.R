#' Find the best bin size for a given dataset
#'
#' Use simulations to find the best bin size among a set of input files. There is no guarantee that the bin size will be the best for your data, since it is only "best" in terms of fewest miscalls for simulated data. However, it can give you a hint what bin size to choose.
#'
#' The function first runs \code{\link{callPeaksUnivariate}} on the given binned.data files. From the estimated parameters it generates simulated data and calls the peaks on this simulated data. Because the data is simulated, the fraction of miscalls can be precisely calculated.
#'
#' @author Aaron Taudt
#' @param files.binned A vector with files that contain \code{\link{binned.data}} in different bin sizes.
#' @param outputfolder Name of the folder where all files will be written to.
#' @param chromosomes A vector of chromosomes to use for the simulation.
#' @param eps Convergence threshold for the Baum-Welch algorithm.
#' @param max.iter The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.
#' @param max.time The maximum running time in seconds for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.
#' @param repetitions Number of repetitions for each simulation.
#' @param plot.progress If TRUE, the plot will be updated each time a simulation has finished. If FALSE, the plot will be returned only at the end.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object with a bar plot of the number of miscalls dependent on the bin size.
scanBinsizes = function(files.binned, outputfolder, chromosomes="chr10", eps=0.01, max.iter=100, max.time=300, repetitions=3, plot.progress=FALSE) {

    ## Create outputfolder if not existent
    if (!file.exists(outputfolder)) dir.create(outputfolder)

    ## Do univariates (real)
    path.uni = file.path(outputfolder, "results_univariate")
    if (!file.exists(path.uni)) dir.create(path.uni)

    for (binfile in files.binned) {
        unifile = file.path(path.uni, paste0("univariate_",basename(binfile)))
        if (!file.exists(unifile)) {
            message()
            messageU(paste0('Univariate for ',basename(binfile),':'))
            binned.data <- get(load(binfile))
            if (!is.null(chromosomes)) {
                binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
            }
            unimodel = callPeaksUnivariate(binned.data, eps=eps, max.iter=max.iter, max.time=max.time)
            save(unimodel, file=unifile)
        }
    }

    message()
    performance = NULL
    for (irep in 1:repetitions) {
        messageU(paste0('Simulation ',irep,':'))

        ## Simulate data from univariates
        path.sim.bin = file.path(outputfolder, paste0("simulated_rep_",irep,"_binned"))
        if (!file.exists(path.sim.bin)) dir.create(path.sim.bin)
        simfile = file.path(path.sim.bin, paste0("simulated.RData"))
        if (!file.exists(simfile)) {
            ptm <- startTimedMessage('  simulating data ...')
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
            stopTimedMessage(ptm)
            # Save with states
            ptm <- startTimedMessage('  saving to file ...')
            save(sim.data.list, file=simfile)
            stopTimedMessage(ptm)
        } else {
            ptm <- startTimedMessage('  loading simulated data ...')
            sim.data.list = get(load(simfile))
            stopTimedMessage(ptm)
        }
        # Get simulated original states
        ptm <- startTimedMessage("  getting simulated original states ...")
        sim.states.list = lapply(lapply(sim.data.list, "[[", 'bins'), function(gr) { return(gr$state) } )
        for (binfile in files.binned) {
            sim.states.list[[basename(binfile)]] = state.labels[sim.states.list[[basename(binfile)]]]
        }
        stopTimedMessage(ptm)

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
                messageU(paste0('Univariate for simulated ',basename(binfile),':'))
                binned.data <- get(load(sim.binfile))
                if (!is.null(chromosomes)) {
                    binned.data <- binned.data[seqnames(binned.data) %in% chromosomes]
                }
                unimodel = callPeaksUnivariate(binned.data, eps=eps, max.iter=max.iter, max.time=max.time)
                save(unimodel, file=sim.unifile)
            } else {
                unimodel <- get(load(sim.unifile))
            }
            # Get univariate predicted states
            ptm <- startTimedMessage("  getting predicted univariate states...")
            uni.states.list[[basename(binfile)]] = unimodel$bins$state
            binsizes[basename(binfile)] <- width(unimodel$bins)[1]
            stopTimedMessage(ptm)

            # Performance
            ptm <- startTimedMessage("  calculating performance...")
            mask = sim.states.list[[basename(binfile)]] != uni.states.list[[basename(binfile)]]
            miscalls = length(which(mask)) / length(sim.states.list[[basename(binfile)]])
            performance = rbind(performance, data.frame(simulation=irep, binsize=binsizes[basename(binfile)], miscalls))
            stopTimedMessage(ptm)

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
