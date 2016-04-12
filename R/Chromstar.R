#' Wrapper function for the \code{\link{chromstaR}} package
#' 
#' This function performs \code{\link[chromstaR:binReads]{binning}}, \code{\link[chromstaR:callPeaksUnivariate]{univariate peak calling}} and \code{\link[chromstaR:callPeaksMultivariate]{multivariate peak calling}} from a list of input files.
#' 
#' @param inputfolder Folder with either BAM or BED files.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param experiment.table A \code{data.frame} or tab-separated text file with the structure of the experiment. See \code{\link{experiment.table}} for an example.
#' @inheritParams binReads
#' @inheritParams callPeaksUnivariate
#' @param mode One of \code{c('condition','mark','full')}. The modes determine how the multivariate part is run. Here is some advice which mode to use:
#' \describe{
#'   \item{\code{mark}}{Each condition is analyzed separately with all marks combined. Choose this mode if you have more than ~7 conditions or you want to have a high sensitivity for detecting combinatorial states. Differences between conditions will be more noisy (more false positives) than in mode \code{'condition'} but combinatorial states are more precise.}
#'   \item{\code{condition}}{Each mark is analyzed separately with all conditions combined. Choose this mode if you are interested in accurate differences. Combinatorial states will be more noisy (more false positives) than in mode \code{'mark'} but differences are more precise.}
#'   \item{\code{full}}{Full analysis of all marks and conditions combined. Best of both, but: Choose this mode only if (number of conditions * number of marks \eqn{\le} 8), otherwise it might be too slow or crash due to memory limitations.}
#' }
#' @param num.states The number of states to use in the multivariate part. If set to \code{NULL}, the maximum number of theoretically possible states is used. CAUTION: This can be very slow or crash if you have too many states. \pkg{\link{chromstaR}} has a built in mechanism to select the best states in case that less states than theoretically possible are specified.
#' @return \code{NULL}
#' @import foreach
#' @import doParallel
#' @importFrom utils read.table
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @export
#' 
#' @examples 
#'## Prepare the file paths. Exchange this with your input and output directories.
#'inputfolder <- system.file("extdata","euratrans", package="chromstaRData")
#'outputfolder <- file.path(tempdir(), 'SHR-example')
#'## Define experiment structure
#'data(experiment_table_SHR)
#'## Define assembly
#'# This is only necessary if you have BED files, BAM files are handled automatically.
#'# For common assemblies you can also specify them as 'hg19' for example.
#'data(rn4_chrominfo)
#'## Run ChromstaR
#'Chromstar(inputfolder, experiment.table=experiment_table_SHR,
#'          outputfolder=outputfolder, numCPU=2, binsize=1000, assembly=rn4_chrominfo,
#'          prefit.on.chr='chr12', mode='mark', eps=1)
#'
Chromstar <- function(inputfolder, experiment.table, outputfolder, configfile=NULL, numCPU=1, binsize=1000, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, prefit.on.chr=NULL, eps=0.01, max.time=NULL, max.iter=5000, read.cutoff.absolute=500, keep.posteriors=FALSE, mode='mark', num.states=128) {
  
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
  params <- list(numCPU=numCPU, binsize=binsize, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, prefit.on.chr=prefit.on.chr, eps=eps, max.time=max.time, max.iter=max.iter, read.cutoff.absolute=read.cutoff.absolute, keep.posteriors=keep.posteriors, mode=mode, num.states=num.states)
  conf <- c(conf, params[setdiff(names(params),names(conf))])
  
  ## Helpers
  binsize <- conf[['binsize']]
  numcpu <- conf[['numCPU']]
  mode <- conf[['mode']]
  
  ## Read in experiment table if necessary ##
	if (is.character(experiment.table)) {
		exp.table <- utils::read.table(experiment.table, header=TRUE, comment.char='#')
		if (!all(colnames(exp.table) == c('file','mark','condition','replicate','pairedEndReads'))) {
			stop("Your 'experiment.table' must be a tab-separated file with column names 'file', 'mark', 'condition', 'replicate' and 'pairedEndReads'.")
		}
		rownames(exp.table) <- exp.table[,1]
	} else if (is.data.frame(experiment.table)) {
		exp.table <- experiment.table
		rownames(exp.table) <- exp.table[,1]
	} else {
		stop("Argument 'experiment.table' must be a data.frame or a tab-separated file.")
	}
  
  ## Check usage of modes
  marks <- unique(as.character(exp.table[,'mark']))
  conditions <- unique(as.character(exp.table[,'condition']))
  if (length(conditions) < 2 & conf[['mode']] == 'condition') {
    stop("Mode 'condition' can only be used if two or more conditions are present.")
  }
  if (length(marks) < 2 & conf[['mode']] == 'mark') {
    stop("Mode 'mark' can only be used if two or more marks are present.")
  }
  
  ## Check if assembly must be present
  files <- file.path(inputfolder, exp.table$file)
  files.clean <- sub('\\.gz$','', files)
  format <- sapply(strsplit(files.clean, '\\.'), function(x) { rev(x)[1] })
  if (any(format=='bed') & is.null(conf[['assembly']])) {
    stop("Please specify an 'assembly' for the BED files.")
  }
  
  ## Read in assembly if necessary
  if (is.character(conf[['assembly']])) {
    if (file.exists(conf[['assembly']])) {
      conf[['assembly']] <- utils::read.table(conf[['assembly']], sep='\t', header=TRUE)
    }
  }
  
  ## Set up the directory structure ##
  binpath <- file.path(outputfolder, 'binned')
  unipath <- file.path(outputfolder, 'univariate')
  plotpath <- file.path(outputfolder, 'plots')
  multipath <- file.path(outputfolder, 'multivariate')
  combipath <- file.path(outputfolder, 'multivariate-combined')
  browserpath <- file.path(outputfolder, 'browserfiles')
  if (!file.exists(outputfolder)) { dir.create(outputfolder) }
  if (!file.exists(plotpath)) { dir.create(plotpath) }
  filenames <- paste0(exp.table$file, "_binsize",format(binsize, scientific=FALSE, trim=TRUE),".RData")
  
  ## Make a copy of the conf file
  writeConfig(conf, configfile=file.path(outputfolder, 'chromstaR.config'))
  
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
  parallel.helper <- function(file) {
    savename <- file.path(binpath, paste0(basename(file),"_binsize",format(binsize, scientific=FALSE, trim=TRUE),".RData"))
    if (!file.exists(savename)) {
      tC <- tryCatch({
        binReads(file=file, assembly=conf[['assembly']], pairedEndReads=exp.table[basename(file),'pairedEndReads'], binsizes=binsize, chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], outputfolder.binned=binpath, save.as.RData=TRUE)
      }, error = function(err) {
        stop(file,'\n',err)
      })
    }
  }
  if (numcpu > 1) {
    ptm <- startTimedMessage("Binning the data ...")
    temp <- foreach (file = files, .packages=c("chromstaR")) %dopar% {
      parallel.helper(file)
    }
    stopTimedMessage(ptm)
  } else {
    temp <- foreach (file = files, .packages=c("chromstaR")) %do% {
      parallel.helper(file)
    }
  }
  
  
  #==============================
  ### Univariate peak calling ###
  #==============================
  if (!file.exists(unipath)) { dir.create(unipath) }
  files <- file.path(binpath, filenames)
  
  parallel.helper <- function(file) {
    savename <- file.path(unipath, basename(file))
    if (!file.exists(savename)) {
      tC <- tryCatch({
        fields <- sapply(exp.table[sub('_binsize.*','',basename(file)),c('mark','condition','replicate')], as.character)
        id <- paste0(fields[1], '-', fields[2], '-rep', fields[3])
        model <- callPeaksUnivariate(file, ID=id, eps=conf[['eps']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], read.cutoff.absolute=conf[['read.cutoff.absolute']], prefit.on.chr=conf[['prefit.on.chr']], keep.posteriors=conf[['keep.posteriors']])
        save(model, file=savename)
      }, error = function(err) {
        stop(file,'\n',err)
      })
    }
  }
  if (numcpu > 1) {
    ptm <- startTimedMessage("Univariate peak calling ...")
    temp <- foreach (file = files, .packages=c("chromstaR")) %dopar% {
      parallel.helper(file)
    }
    stopTimedMessage(ptm)
  } else {
    temp <- foreach (file = files, .packages=c("chromstaR")) %do% {
      parallel.helper(file)
    }
  }
  
  #-----------------------
  ## Plot distributions ##
  #-----------------------
  # ptm <- startTimedMessage("Plotting univariate distributions ...")
  # files <- file.path(unipath, filenames)
  # 
  # savename <- file.path(plotpath, 'univariate-distributions.pdf')
  # grDevices::pdf(savename, width=7, height=5)
  # for (file in files) {
  #   print(graphics::plot(file))
  # }
  # d <- grDevices::dev.off()
  # stopTimedMessage(ptm)
  
  
  #================================
  ### Multivariate peak calling ###
  #================================
  if (!file.exists(multipath)) { dir.create(multipath) }
  if (!file.exists(browserpath)) { dir.create(browserpath) }

	## Plot helper ##
	plothelper <- function(savename, multimodel) {
	  char.per.cm <- 10
	  legend.cm <- 3
		savename1 <- paste0(savename, '_correlation.pdf')
		ggplt <- graphics::plot(multimodel, type='correlation')
		width <- length(multimodel$IDs) + max(sapply(multimodel$IDs, nchar)) / char.per.cm + legend.cm
		height <- length(multimodel$IDs) + max(sapply(multimodel$IDs, nchar)) / char.per.cm
		ggsave(savename1, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')

		savename2 <- paste0(savename, '_transitionMatrix.pdf')
		ggplt <- graphics::plot(multimodel, type='transitionMatrix')
		width <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + legend.cm
		height <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + 1
		ggsave(savename2, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
	}
  
  ## Run multivariate depending on mode
  multimodels <- list()
	#--------------------
  if (mode == 'full') {
    savename <- file.path(multipath, paste0('multivariate_mode-', mode, '.RData'))
    if (!file.exists(savename)) {
      files <- file.path(unipath, filenames)
      states <- stateBrewer(exp.table, mode=mode)
      multimodel <- callPeaksMultivariate(files, use.states=states, num.states=conf[['num.states']], eps=conf[['eps']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']])
      save(multimodel, file=savename)
    } else {
      multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
    }
		multimodels[[1]] <- multimodel
    ## Export browser files
    savename <- file.path(browserpath, paste0('multivariate_mode-', mode))
    if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
      trackname <- paste0('combinations, mode-', mode)
      exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
    }
    if (!file.exists(paste0(savename, '_counts.wig.gz'))) {
      exportMultivariate(multimodel, filename=savename, what='counts')
    }
		## Plot transition and correlations
		savename <- file.path(plotpath, paste0('multivariate_mode-', mode))
		plothelper(savename, multimodel)
    
	#---------------------------
  } else if (mode == 'mark') {
    for (condition in conditions) {
      savename <- file.path(multipath, paste0('multivariate_mode-', mode, '_condition-', condition, '.RData'))
      if (!file.exists(savename)) {
        mask <- exp.table[,'condition'] == condition
        files <- file.path(unipath, filenames)[mask]
        states <- stateBrewer(exp.table[mask,], mode=mode)
        multimodel <- callPeaksMultivariate(files, use.states=states, num.states=conf[['num.states']], eps=conf[['eps']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']])
        save(multimodel, file=savename)
      } else {
				multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
      }
			multimodels[[as.character(condition)]] <- multimodel
      ## Export browser files
      savename <- file.path(browserpath, paste0('multivariate_mode-', mode, '_condition-', condition))
      if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        trackname <- paste0('combinations, mode-', mode)
        exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
      }
      if (!file.exists(paste0(savename, '_counts.wig.gz'))) {
        exportMultivariate(multimodel, filename=savename, what='counts')
      }
			## Plot transition and correlations
			savename <- file.path(plotpath, paste0('multivariate_mode-', mode, '_condition-', condition))
			plothelper(savename, multimodel)
    }
    
	#--------------------------------
  } else if (mode == 'condition') {
    for (mark in marks) {
      savename <- file.path(multipath, paste0('multivariate_mode-', mode, '_mark-', mark, '.RData'))
      if (!file.exists(savename)) {
        mask <- exp.table[,'mark'] == mark
        files <- file.path(unipath, filenames)[mask]
        states <- stateBrewer(exp.table[mask,], mode=mode)
        multimodel <- callPeaksMultivariate(files, use.states=states, num.states=conf[['num.states']], eps=conf[['eps']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']])
        save(multimodel, file=savename)
      } else {
				multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
      }
			multimodels[[as.character(mark)]] <- multimodel
      ## Export browser files
      savename <- file.path(browserpath, paste0('multivariate_mode-', mode, '_mark-', mark))
      if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        trackname <- paste0('combinations, mode-', mode)
        exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
      }
      if (!file.exists(paste0(savename, '_counts.wig.gz'))) {
        exportMultivariate(multimodel, filename=savename, what='counts')
      }
			## Plot transition and correlations
			savename <- file.path(plotpath, paste0('multivariate_mode-', mode, '_mark-', mark))
			plothelper(savename, multimodel)
    }
  }

  
  #================================
  ## Combine multiple conditions ##
  #================================
  if ((mode=='mark' & length(conditions)>=2) | (mode=='condition' & length(marks)>=2)) {
    
    if (!file.exists(combipath)) { dir.create(combipath) }
    savename <- file.path(combipath, paste0('combined_mode-', mode, '.RData'))
    if (!file.exists(savename)) {
      combinedModel <- combineMultivariates(multimodels, mode=mode, conditions=conditions)
      save(combinedModel, file=savename)
    } else {
			combinedModel <- loadHmmsFromFiles(savename, check.class=class.combined.multivariate.hmm)[[1]]
    }
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    if (!file.exists(browserpath)) { dir.create(browserpath) }
    savename <- file.path(browserpath, paste0('combined_mode-', mode))
    if (!file.exists(paste0(savename,'.bed.gz'))) {
      trackname <- paste0('combinations, mode-', mode)
      exportCombinedMultivariate(combinedModel, filename=savename, trackname=trackname)
    }
  }
  
  total.time <- proc.time() - total.time
  message("==> Total time spent: ", round(total.time[3]), "s <==")
  
  
}
