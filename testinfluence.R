load_all()

savename <- '~/Bioconductor/chromstaR/influencetest.RData'
if (!file.exists(savename)) {
  # Get example BAM files for 2 different marks in hypertensive rat
  file.path <- system.file("extdata","euratrans", package='chromstaRData')
  files <- list.files(file.path, full.names=TRUE, pattern='SHR.*bam$')[c(1:2,6)]
  # Construct experiment structure
  exp <- data.frame(file=files, mark=c("H3K27me3","H3K27me3","H3K4me3"),
                    condition=rep("SHR",3), replicate=c(1:2,1), pairedEndReads=FALSE,
                    controlFiles=NA)
  states <- stateBrewer(exp, mode='combinatorial')
  # Bin the data
  data(rn4_chrominfo)
  binned.data <- list()
  for (file in files) {
    binned.data[[basename(file)]] <- binReads(file, binsizes=1000, stepsizes=500,
                                              experiment.table=exp,
                                              assembly=rn4_chrominfo, chromosomes='chr12')
  }
  # Obtain the univariate fits
  models <- list()
  for (i1 in 1:length(binned.data)) {
    models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, eps=1)
  }
  save(models, states, file = savename)
} else {
  load(savename)
}
  
# Call multivariate peaks
load_all(); multimodel <- callPeaksInfluence(models, use.states=states, eps=1, max.time=60, verbosity=1, keep.densities = TRUE, max.iter=NULL)
