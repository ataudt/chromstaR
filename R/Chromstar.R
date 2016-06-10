#' Wrapper function for the \code{\link{chromstaR}} package
#' 
#' This function performs \code{\link[chromstaR:binReads]{binning}}, \code{\link[chromstaR:callPeaksUnivariate]{univariate peak calling}} and \code{\link[chromstaR:callPeaksMultivariate]{multivariate peak calling}} from a list of input files.
#' 
#' @param inputfolder Folder with either BAM or BED files.
#' @param experiment.table A \code{data.frame} or tab-separated text file with the structure of the experiment. See \code{\link{experiment.table}} for an example.
#' @param outputfolder Folder where the results and intermediate files will be written to.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU Number of threads to use for the analysis. Beware that more CPUs also means more memory is needed. If you experience crashes of R with higher numbers of this parameter, leave it at \code{numCPU=1}.
#' @param binsize An integer specifying the bin size that is used for the analysis.
#' @inheritParams readBedFileAsGRanges
#' @inheritParams callPeaksUnivariate
#' @param mode One of \code{c('condition','mark','full')}. The modes determine how the multivariate part is run. Here is some advice which mode to use:
#' \describe{
#'   \item{\code{mark}}{Each condition is analyzed separately with all marks combined. Choose this mode if you have more than ~7 conditions or you want to have a high sensitivity for detecting combinatorial states. Differences between conditions will be more noisy (more false positives) than in mode \code{'condition'} but combinatorial states are more precise.}
#'   \item{\code{condition}}{Each mark is analyzed separately with all conditions combined. Choose this mode if you are interested in accurate differences. Combinatorial states will be more noisy (more false positives) than in mode \code{'mark'} but differences are more precise.}
#'   \item{\code{full}}{Full analysis of all marks and conditions combined. Best of both, but: Choose this mode only if (number of conditions * number of marks \eqn{\le} 8), otherwise it might be too slow or crash due to memory limitations.}
#' }
#' @param max.states The maximum number of states to use in the multivariate part. If set to \code{NULL}, the maximum number of theoretically possible states is used. CAUTION: This can be very slow or crash if you have too many states. \pkg{\link{chromstaR}} has a built in mechanism to select the best states in case that less states than theoretically possible are specified.
#' @param per.chrom If set to \code{TRUE} chromosomes will be treated separately in the multivariate part. This tremendously speeds up the calculation but results might be noisier as compared to \code{per.chrom=FALSE}, where all chromosomes are concatenated for the HMM.
#' @param eps.univariate Convergence threshold for the univariate Baum-Welch algorithm.
#' @param eps.multivariate Convergence threshold for the multivariate Baum-Welch algorithm.
#' @return \code{NULL}
#' @import foreach
#' @import doParallel
#' @importFrom utils read.table write.table
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
#'          prefit.on.chr='chr12', mode='mark', eps.univariate=1, eps.multivariate=1)
#'
Chromstar <- function(inputfolder, experiment.table, outputfolder, configfile=NULL, numCPU=1, binsize=1000, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, prefit.on.chr=NULL, eps.univariate=0.1, max.time=NULL, max.iter=5000, read.cutoff.absolute=500, keep.posteriors=TRUE, mode='mark', max.states=128, per.chrom=TRUE, eps.multivariate=0.01) {
  
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
    params <- list(numCPU=numCPU, binsize=binsize, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, prefit.on.chr=prefit.on.chr, eps.univariate=eps.univariate, max.time=max.time, max.iter=max.iter, read.cutoff.absolute=read.cutoff.absolute, keep.posteriors=keep.posteriors, mode=mode, max.states=max.states, per.chrom=per.chrom, eps.multivariate=eps.multivariate)
    conf <- c(conf, params[setdiff(names(params),names(conf))])
    
    ## Helpers
    binsize <- conf[['binsize']]
    binsize.string <- format(binsize, scientific=FALSE, trim=TRUE)
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
    if (!mode %in% c('mark','condition','full')) {
        stop("Unknown mode '", mode, "'.")
    }
    marks <- setdiff(unique(as.character(exp.table[,'mark'])), 'input')
    conditions <- unique(as.character(exp.table[,'condition']))
    if (length(conditions) < 2 & conf[['mode']] == 'condition') {
        stop("Mode 'condition' can only be used if two or more conditions are present.")
    }
    if (length(marks) < 2 & conf[['mode']] == 'mark') {
        stop("Mode 'mark' can only be used if two or more marks are present.")
    }
    
    ## Check if assembly must be present
    datafiles <- file.path(inputfolder, basename(as.character(exp.table$file)))
    datafiles.clean <- sub('\\.gz$','', datafiles)
    format <- sapply(strsplit(datafiles.clean, '\\.'), function(x) { rev(x)[1] })
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
    reppath <- file.path(outputfolder, 'replicates')
    plotpath <- file.path(outputfolder, 'PLOTS')
    uniplotpath <- file.path(plotpath, 'univariate-distributions')
    multipath <- file.path(outputfolder, 'multivariate')
    combipath <- file.path(outputfolder, 'combined')
    browserpath <- file.path(outputfolder, 'BROWSERFILES')
    if (!file.exists(outputfolder)) { dir.create(outputfolder) }
    if (!file.exists(plotpath)) { dir.create(plotpath) }
    if (!file.exists(browserpath)) { dir.create(browserpath) }
    
    ## Write README ##
    savename <- file.path(outputfolder, 'README.txt')
    cat("", file=savename)
    cat("This folder contains the following files:\n", file=savename, append=TRUE)
    cat("-----------------------------------------\n", file=savename, append=TRUE)
    cat("- chrominfo.tsv: A tab-separated file with chromosome lengths.\n", file=savename, append=TRUE)
    cat("- chromstaR.config: A text file with all the parameters that were used to run Chromstar().\n", append=TRUE, file=savename)
    cat("- experiment_table.tsv: A tab-separated file of your experiment setup.\n", file=savename, append=TRUE)
    
    cat("\n", file=savename, append=TRUE)
    cat("This folder contains the following folders:\n", file=savename, append=TRUE)
    cat("-------------------------------------------\n", file=savename, append=TRUE)
    cat("- binned: RData files with the results of the binnig step. Contains GRanges objects with binned genomic coordinates and read counts.\n", file=savename, append=TRUE)
    cat("- BROWSERFILES: Bed files for upload to the UCSC genome browser. It contains files with combinatorial states (*_combinations.bed.gz), underlying peak calls (*_peaks.bed.gz), and read counts (*_counts.wig.gz). !!Always check the *_peaks.bed.gz files if you are satisfied with the peak calls. If not, there are ways to make the calls stricter (see section FAQ of the vignette).\n", file=savename, append=TRUE)
    cat("- -->combined<--: RData files with the combined results of the uni- and multivariate peak calling steps. This is what you want to use for downstream analyses. Contains combinedMultiHMM objects.\n", file=savename, append=TRUE)
    cat("    - combined_mode-separate.RData: Simple combination of univariate peak calls (replicates considered) without multivariate analysis.\n", file=savename, append=TRUE)
    cat("    - combined_mode-mark.RData: Combination of multivariate results for mode='mark'.\n", file=savename, append=TRUE)
    cat("    - combined_mode-condition.RData: Combination of multivariate results for mode='condition'.\n", file=savename, append=TRUE)
    cat("    - combined_mode-full.RData: Combination of multivariate results for mode='full'.\n", file=savename, append=TRUE)
    cat("- multivariate: RData files with the results of the multivariate peak calling step. Contains multiHMM objects.\n", file=savename, append=TRUE)
    cat("- PLOTS: Several plots that are produced by default. Please check the plots in subfolder \"univariate-distributions\" for irregularities (see section \"Univariate Analysis\" of the vignette).\n", file=savename, append=TRUE)
    cat("- replicates: RData files with the result of the replicate peak calling step. Contains multiHMM objects.\n", file=savename, append=TRUE)
    cat("- univariate: RData files with the result of the univariate peak calling step. Contains uniHMM objects.\n", file=savename, append=TRUE)
    
    ## Make a copy of the conf file
    writeConfig(conf, configfile=file.path(outputfolder, 'chromstaR.config'))

    ## Write the experiment table to file
    utils::write.table(exp.table, file=file.path(outputfolder, 'experiment_table.tsv'), col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')

    ## Prepare IDs and filenames
    IDs <- paste0(exp.table$mark, '-', exp.table$condition, '-rep', exp.table$replicate)
    names(IDs) <- basename(datafiles)
    filenames <- paste0(IDs, '_binsize', binsize.string, '.RData')
    names(filenames) <- basename(datafiles)

    
  
    #==============
    ### Binning ###
    #==============
    if (!file.exists(binpath)) { dir.create(binpath) }
    parallel.helper <- function(file) {
        savename <- file.path(binpath, filenames[basename(file)])
        if (!file.exists(savename)) {
            tC <- tryCatch({
                binReads(file=file, experiment.table=exp.table, assembly=conf[['assembly']], pairedEndReads=exp.table[basename(file),'pairedEndReads'], binsizes=binsize, chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], outputfolder.binned=binpath, save.as.RData=TRUE)
            }, error = function(err) {
                stop(file,'\n',err)
            })
        }
    }
    for (file in datafiles) {
        parallel.helper(file)
    }
    
    ## Delete input rows from exp.table
    exp.table <- exp.table[exp.table$mark != 'input',]
    IDs <- IDs[!grepl('input', IDs)]
    filenames <- filenames[!grepl('input', basename(filenames))]
  
  
    #==============================
    ### Univariate peak calling ###
    #==============================
    ## Parallelization ##
    if (numcpu > 1) {
        ptm <- startTimedMessage("Setting up parallel execution with ", numcpu, " threads ...")
        cl <- parallel::makeCluster(numcpu)
        doParallel::registerDoParallel(cl)
        on.exit(
            if (conf[['numCPU']] > 1) {
                parallel::stopCluster(cl)
            }
        )
        stopTimedMessage(ptm)
    }
  
    if (!file.exists(unipath)) { dir.create(unipath) }
    if (!file.exists(uniplotpath)) { dir.create(uniplotpath) }
    files <- file.path(binpath, filenames)
    inputfiles <- files[grep('input', basename(files))]
    files <- files[grep('input', basename(files), invert=TRUE)]
    
    parallel.helper <- function(file, inputfiles) {
        ## Peak calling
        savename <- file.path(unipath, basename(file))
        if (!file.exists(savename)) {
            tC <- tryCatch({
                mark <- strsplit(basename(file), '-')[[1]][1]
                condition <- strsplit(basename(file), '-')[[1]][2]
                input.files <- grep(paste0('input-',condition,'-'), inputfiles, value=TRUE)
                if (length(input.files) == 0) {
                    input.files <- NULL
                }
                model <- callPeaksUnivariate(file, input.data=input.files, eps=conf[['eps.univariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], read.cutoff.absolute=conf[['read.cutoff.absolute']], prefit.on.chr=conf[['prefit.on.chr']], keep.posteriors=FALSE)
                ptm <- startTimedMessage("Saving to file ", savename, " ...")
                save(model, file=savename)
                stopTimedMessage(ptm)
            }, error = function(err) {
                stop(file,'\n',err)
            })
        }
        ## Plot distributions
        savename.pdf <- file.path(uniplotpath, paste0(basename(file),'.pdf'))
        if (!file.exists(savename.pdf)) {
            ptm <- startTimedMessage("Plotting to file ", savename.pdf, " ...")
            ggplt <- plotHistogram(savename)
            ggsave(filename=savename.pdf, plot=ggplt, width=7, height=5)
            stopTimedMessage(ptm)
        }
    }
    if (numcpu > 1) {
        ptm <- startTimedMessage("Univariate peak calling ...")
        temp <- foreach (file = files, .packages=c("chromstaR")) %dopar% {
            parallel.helper(file, inputfiles)
        }
        stopTimedMessage(ptm)
    } else {
        for (file in files) {
            parallel.helper(file, inputfiles)
        }
    }
  
  
    #================================
    ### Replicates                ###
    #================================
    messageU("Analyzing replicates")
    if (!file.exists(reppath)) { dir.create(reppath) }
    
    files <- file.path(unipath, filenames)
    names(files) <- paste0(exp.table$mark, '-', exp.table$condition)
    markconditions <- unique(names(files))
    for (markcond in markconditions) {
        savename <- file.path(reppath, paste0(markcond, '_binsize', binsize, '.RData'))
        if (!file.exists(savename)) {
            mask <- names(files) == markcond
            repfiles <- files[mask]
            states <- stateBrewer(exp.table[mask,], mode='mark')
            repmodel <- callPeaksMultivariate(repfiles, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']])
            ptm <- startTimedMessage("Saving to file ", savename, " ...")
            save(repmodel, file=savename)
            stopTimedMessage(ptm)
        } else {
            repmodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
        }
    }
    
    message("Combining replicate HMMs")
    modename <- 'separate'
    repfiles <- list.files(reppath, full.names=TRUE)
    ## Get the order correct
    names(repfiles) <- gsub('_binsize.*', '', basename(repfiles))
    ordering <- unique(gsub('-rep.*', '', IDs))
    repfiles <- repfiles[ordering]
    if (!file.exists(combipath)) { dir.create(combipath) }
    savename <- file.path(combipath, paste0('combined_mode-', modename, '.RData'))
    if (!file.exists(savename)) {
        combined.model <- combineMultivariates(repfiles, mode='replicate')
        ptm <- startTimedMessage("Saving to file ", savename, " ...")
        save(combined.model, file=savename)
        stopTimedMessage(ptm)
    } else {
        combined.model <- loadHmmsFromFiles(savename, check.class=class.combined.multivariate.hmm)[[1]]
    }
    ## Plot correlations
    ptm <- startTimedMessage("Plotting read count correlation ...")
    char.per.cm <- 10
    legend.cm <- 3
    savename <- file.path(plotpath, 'read-count-correlation.pdf')
    ggplt <- heatmapCountCorrelation(combined.model, cluster=FALSE)
    width <- length(combined.model$info$ID) + max(sapply(combined.model$info$ID, nchar)) / char.per.cm + legend.cm
    height <- length(combined.model$info$ID) + max(sapply(combined.model$info$ID, nchar)) / char.per.cm
    ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
    savename <- file.path(plotpath, 'read-count-correlation-clustered.pdf')
    ggplt <- heatmapCountCorrelation(combined.model, cluster=TRUE)
    width <- length(combined.model$info$ID) + max(sapply(combined.model$info$ID, nchar)) / char.per.cm + legend.cm
    height <- length(combined.model$info$ID) + max(sapply(combined.model$info$ID, nchar)) / char.per.cm
    ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    if (!file.exists(browserpath)) { dir.create(browserpath) }
    savename <- file.path(browserpath, paste0('combined_mode-', modename))
    trackname <- paste0('mode-', modename)
    if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='combinations')
    }
    if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='peaks')
    }
    

    #================================
    ### Multivariate peak calling ###
    #================================
    messageU("Calling multivariate peaks")
    message("mode = ", mode)
    if (!file.exists(multipath)) { dir.create(multipath) }
    if (!file.exists(browserpath)) { dir.create(browserpath) }

    ## Plot helper ##
    plothelper <- function(savename, multimodel) {
        ptm <- startTimedMessage("Making plots ...")
        char.per.cm <- 10
        legend.cm <- 3

        savename2 <- paste0(savename, '_transitionMatrix.pdf')
        ggplt <- suppressMessages( heatmapTransitionProbs(multimodel) )
        width <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + legend.cm
        height <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + 1
        ggsave(savename2, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
        stopTimedMessage(ptm)
    }
  
    ## Run multivariate depending on mode
    #--------------------
    if (mode == 'full') {
        savename <- file.path(multipath, paste0('multivariate_mode-', mode, '.RData'))
        if (!file.exists(savename)) {
            files <- file.path(unipath, filenames)
            states <- stateBrewer(exp.table, mode=mode)
            multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']])
            ptm <- startTimedMessage("Saving to file ", savename, " ...")
            save(multimodel, file=savename)
            stopTimedMessage(ptm)
        } else {
            multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
        }
        ## Export browser files
        savename <- file.path(browserpath, paste0('multivariate_mode-', mode))
        # if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        #     trackname <- paste0('combinations, mode-', mode)
        #     exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
        # }
        # if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
        #     exportMultivariate(multimodel, filename=savename, what='peaks')
        # }
        if (!file.exists(paste0(savename, '_counts.wig.gz'))) {
            exportMultivariate(multimodel, filename=savename, what='counts')
        }
        ## Plot transition and correlations
        savename <- file.path(plotpath, paste0('multivariate_mode-', mode))
        plothelper(savename, multimodel)
    
    #---------------------------
    } else if (mode == 'mark') {
        for (condition in conditions) {
            messageU("condition = ", condition, underline='-')
            savename <- file.path(multipath, paste0('multivariate_mode-', mode, '_condition-', condition, '.RData'))
            if (!file.exists(savename)) {
                mask <- exp.table[,'condition'] == condition
                files <- file.path(unipath, filenames)[mask]
                states <- stateBrewer(exp.table[mask,], mode=mode)
                multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']])
                ptm <- startTimedMessage("Saving to file ", savename, " ...")
                save(multimodel, file=savename)
                stopTimedMessage(ptm)
            } else {
                multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
            }
            ## Export browser files
            savename <- file.path(browserpath, paste0('multivariate_mode-', mode, '_condition-', condition))
            # if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
            #     trackname <- paste0('combinations, mode-', mode)
            #     exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
            # }
            # if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
            #     exportMultivariate(multimodel, filename=savename, what='peaks')
            # }
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
            messageU("mark = ", mark, underline='-')
            savename <- file.path(multipath, paste0('multivariate_mode-', mode, '_mark-', mark, '.RData'))
            if (!file.exists(savename)) {
                mask <- exp.table[,'mark'] == mark
                files <- file.path(unipath, filenames)[mask]
                states <- stateBrewer(exp.table[mask,], mode=mode)
                multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']])
                ptm <- startTimedMessage("Saving to file ", savename, " ...")
                save(multimodel, file=savename)
                stopTimedMessage(ptm)
            } else {
                multimodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
            }
            ## Export browser files
            savename <- file.path(browserpath, paste0('multivariate_mode-', mode, '_mark-', mark))
            # if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
            #     trackname <- paste0('combinations, mode-', mode)
            #     exportMultivariate(multimodel, filename=savename, what='combinations', trackname=trackname)
            # }
            # if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
            #     exportMultivariate(multimodel, filename=savename, what='peaks')
            # }
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
    messageU("Combining multivariate HMMs")
    multifiles <- list.files(multipath, full.names=TRUE, pattern=paste0('mode-',mode))
    if (mode=='condition') {
        names(multifiles) <- marks
    }
    if (mode=='mark') {
        names(multifiles) <- conditions
    }

    if (!file.exists(combipath)) { dir.create(combipath) }
    savename <- file.path(combipath, paste0('combined_mode-', mode, '.RData'))
    if (!file.exists(savename)) {
        combined.model <- combineMultivariates(multifiles, mode=mode)
        ptm <- startTimedMessage("Saving to file ", savename, " ...")
        save(combined.model, file=savename)
        stopTimedMessage(ptm)
    } else {
        combined.model <- loadHmmsFromFiles(savename, check.class=class.combined.multivariate.hmm)[[1]]
    }
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    if (!file.exists(browserpath)) { dir.create(browserpath) }
    savename <- file.path(browserpath, paste0('combined_mode-', mode))
    trackname <- paste0('mode-', mode)
    if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='combinations')
    }
    if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='peaks')
    }
  
    total.time <- proc.time() - total.time
    message("==> Total time spent: ", round(total.time[3]), "s <==")
  
  
}
