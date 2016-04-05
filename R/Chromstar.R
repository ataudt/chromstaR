#' Wrapper function for the \code{\link{chromstaR}} package
#' 
#' This function performs \code{\link[chromstaR:binReads]{binning}}, \code{\link[chromstaR:callPeaksUnivariate]{univariate peak calling}} and \code{\link[chromstaR:callPeaksMultivariate]{multivariate peak calling}} from a list of input files.
#' 
#' @param inputfolder Folder with either BAM or BED files.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param experiment.table A \code{data.frame} or tab-separated text file with the structure of the experiment. See \code{\link{experiment.table}} for an example.
#' @inheritParams binReads
#' @inheritParams callPeaksUnivariate
Chromstar <- function(inputfolder, experiment.table, outputfolder, configfile=NULL, numCPU=1, binsize=1000, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, prefit.on.chr=NULL, eps=0.01, max.time=NULL, max.iter=5000, read.cutoff.absolute=500, keep.posteriors=FALSE) {
  
  #=======================
  ### Helper functions ###
  #=======================
  as.object <- function(x) {
    return(eval(parse(text=x)))
  }
  
  #========================
  ### General variables ###
  #========================
  conf <- NULL
  if (is.character(configfile)) {
    ## Read config file ##
    errstring <- tryCatch({
      conf <- readConfig(configfile)
      errstring <- ''
    }, error = function(err) {
      errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
      stop(errstring)
    }
  }
  total.time <- proc.time()
  
  ## Put options into list and merge with conf
  params <- list(numCPU=numCPU, binsize=binsize, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, prefit.on.chr=prefit.on.chr, eps=eps, max.time=max.time, max.iter=max.iter, read.cutoff.absolute=read.cutoff.absolute, keep.posteriors=keep.posteriors)
  conf <- c(conf, params[setdiff(names(params),names(conf))])
  
  ## Helpers
  binsize <- conf[['binsize']]
  numcpu <- conf[['numCPU']]
  
  ## Set up the directory structure ##
  binpath <- file.path(outputfolder, 'binned')
  unipath <- file.path(outputfolder, 'univariate')
  multipath <- file.path(outputfolder, 'multivariate')
  
  ## Make a copy of the conf file
  writeConfig(conf, configfile=file.path(outputfolder, 'AneuFinder.config'))
  
  ## Read experiment table ##
  exp.table <- read.table(conf[['experiment.table']], header=TRUE, comment.char='#')
  
  ## Parallelization ##
  if (numcpu > 1) {
    ptm <- startTimedMessage("Setting up parallel execution with ", numcpu, " CPUs ...")
    cl <- parallel::makeCluster(numcpu)
    doParallel::registerDoParallel(cl)
    on.exit(
      if (conf[['numCPU']] > 1) {
        parallel::stopCluster(cl)
      }
    )
    stopTimedMessage(ptm)
  }
  
  
  #==============
  ### Binning ###
  #==============
  if (!file.exists(binpath)) { dir.create(binpath) }
  files <- file.path(inputfolder, exp.table$file)
  
  ### Binning ###
  parallel.helper <- function(file) {
    existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value=TRUE)
    existing.binsizes <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 2)))
    existing.rpbin <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 3)))
    binsizes.todo <- setdiff(binsizes, existing.binsizes)
    rpbin.todo <- setdiff(reads.per.bins, existing.rpbin)
    if (length(c(binsizes.todo,rpbin.todo)) > 0) {
      tC <- tryCatch({
        binReads(file=file, format=conf[['format']], assembly=chrom.lengths.df, pairedEndReads=conf[['pairedEndReads']], binsizes=NULL, variable.width.reference=NULL, reads.per.bin=rpbin.todo, bins=bins[as.character(binsizes.todo)], stepsize=conf[['stepsize']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']], outputfolder.binned=binpath.uncorrected, save.as.RData=TRUE, reads.store=TRUE, outputfolder.reads=readspath)
      }, error = function(err) {
        stop(file,'\n',err)
      })
    }
  }
  if (numcpu > 1) {
    ptm <- startTimedMessage("Binning the data ...")
    temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
      parallel.helper(file)
    }
    stopTimedMessage(ptm)
  } else {
    temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
      parallel.helper(file)
    }
  }
  
  
  
  
  
  
  
}