#' Wrapper function for the \pkg{\link{chromstaR}} package
#' 
#' This function performs \code{\link[chromstaR:binReads]{binning}}, \code{\link[chromstaR:callPeaksUnivariate]{univariate peak calling}} and \code{\link[chromstaR:callPeaksMultivariate]{multivariate peak calling}} from a list of input files.
#' 
#' @param inputfolder Folder with either BAM or BED-6 (see \code{\link{readBedFileAsGRanges}} files.
#' @param experiment.table A \code{data.frame} or tab-separated text file with the structure of the experiment. See \code{\link{experiment.table}} for an example.
#' @param outputfolder Folder where the results and intermediate files will be written to.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU Number of threads to use for the analysis. Beware that more CPUs also means more memory is needed. If you experience crashes of R with higher numbers of this parameter, leave it at \code{numCPU=1}.
#' @param binsize An integer specifying the bin size that is used for the analysis.
#' @param stepsize An integer specifying the step size for analysis.
#' @param assembly A \code{data.frame} or tab-separated file with columns 'chromosome' and 'length'. Alternatively a character specifying the assembly, see \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} for available assemblies. Specifying an assembly is only necessary when importing BED files. BAM files are handled automatically.
#' @inheritParams readBedFileAsGRanges
#' @inheritParams callPeaksUnivariate
#' @param mode One of \code{c('differential','combinatorial','full')}. The modes determine how the multivariate part is run. Here is some advice which mode to use:
#' \describe{
#'   \item{\code{combinatorial}}{Each condition is analyzed separately with all marks combined. Choose this mode if you have more than ~7 conditions or you want to have a high sensitivity for detecting combinatorial states. Differences between conditions will be more noisy (more false positives) than in mode \code{'differential'} but combinatorial states are more precise.}
#'   \item{\code{differential}}{Each mark is analyzed separately with all conditions combined. Choose this mode if you are interested in accurate differences. Combinatorial states will be more noisy (more false positives) than in mode \code{'combinatorial'} but differences are more precise.}
#'   \item{\code{full}}{Full analysis of all marks and conditions combined. Best of both, but: Choose this mode only if (number of conditions * number of marks \eqn{\le} 8), otherwise it might be too slow or crash due to memory limitations.}
#'   \item{\code{separate}}{Only replicates are analyzed multivariately. Combinatorial states are constructed by a simple post-hoc combination of peak calls.}
#' }
#' @param max.states The maximum number of states to use in the multivariate part. If set to \code{NULL}, the maximum number of theoretically possible states is used. CAUTION: This can be very slow or crash if you have too many states. \pkg{\link{chromstaR}} has a built in mechanism to select the best states in case that less states than theoretically possible are specified.
#' @param per.chrom If set to \code{TRUE} chromosomes will be treated separately in the multivariate part. This tremendously speeds up the calculation but results might be noisier as compared to \code{per.chrom=FALSE}, where all chromosomes are concatenated for the HMM.
#' @param eps.univariate Convergence threshold for the univariate Baum-Welch algorithm.
#' @param eps.multivariate Convergence threshold for the multivariate Baum-Welch algorithm.
#' @param exclusive.table A \code{data.frame} or tab-separated file with columns 'mark' and 'group'. Histone marks with the same group will be treated as mutually exclusive.
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
#'          outputfolder=outputfolder, numCPU=4, binsize=1000, assembly=rn4_chrominfo,
#'          prefit.on.chr='chr12', chromosomes='chr12', mode='combinatorial', eps.univariate=1,
#'          eps.multivariate=1)
#'
Chromstar <- function(inputfolder, experiment.table, outputfolder, configfile=NULL, numCPU=1, binsize=1000, stepsize=binsize/2, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, prefit.on.chr=NULL, eps.univariate=0.1, max.time=NULL, max.iter=5000, read.cutoff.absolute=500, keep.posteriors=TRUE, mode='differential', max.states=128, per.chrom=TRUE, eps.multivariate=0.01, exclusive.table=NULL) {
  
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
    params <- list(numCPU=numCPU, binsize=binsize, stepsize=stepsize, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, prefit.on.chr=prefit.on.chr, eps.univariate=eps.univariate, max.time=max.time, max.iter=max.iter, read.cutoff.absolute=read.cutoff.absolute, keep.posteriors=keep.posteriors, mode=mode, max.states=max.states, per.chrom=per.chrom, eps.multivariate=eps.multivariate, exclusive.table=exclusive.table)
    conf <- c(conf, params[setdiff(names(params),names(conf))])
    
    ## Helpers
    binsize <- conf[['binsize']]
    stepsize <- conf[['stepsize']]
    binsize.string <- format(binsize, scientific=FALSE, trim=TRUE)
    stepsize.string <- format(stepsize, scientific=FALSE, trim=TRUE)
    numcpu <- conf[['numCPU']]
    mode <- conf[['mode']]
  
    ## Read in experiment table if necessary ##
    if (is.character(experiment.table)) {
        exp.table <- utils::read.table(experiment.table, header=TRUE, comment.char='#', stringsAsFactors=FALSE)
    } else if (is.data.frame(experiment.table)) {
        exp.table <- experiment.table
    }
    check.experiment.table(exp.table)
    
    ## Prepare IDs and filenames
    IDs <- paste0(exp.table$mark, '-', exp.table$condition, '-rep', exp.table$replicate)
    datafiles <- file.path(inputfolder, basename(as.character(exp.table$file)))
    names(datafiles) <- IDs
    rownames(exp.table) <- IDs
    filenames <- paste0(IDs, '_', as.character(exp.table$file), '_binsize', binsize.string, '_stepsize', stepsize.string, '.RData')
    names(filenames) <- IDs
    filenames.stepsize <- paste0(IDs, '_', as.character(exp.table$file), '_binsize', stepsize.string, '_stepsize', stepsize.string, '.RData') # filenames for saving counts at stepsize
    names(filenames.stepsize) <- IDs
    ## Inputfiles
    inputfiles <- unique(unlist(strsplit(as.character(exp.table$controlFiles), '\\|')))
    inputfiles <- inputfiles[!is.na(inputfiles)]
    inputfiles <- file.path(inputfolder, basename(as.character(inputfiles)))
    inputfilenames <- paste0(basename(inputfiles), '_binsize', binsize.string, '_stepsize', stepsize.string, '.RData')
    names(inputfilenames) <- basename(inputfiles)
    inputfilenames.stepsize <- paste0(basename(inputfiles), '_binsize', stepsize.string, '_stepsize', stepsize.string, '.RData')
    names(inputfilenames.stepsize) <- basename(inputfiles)
    
    ## Check usage of modes
    if (!mode %in% c('separate','combinatorial','differential','full')) {
        stop("Unknown mode '", mode, "'.")
    }
    marks <- setdiff(unique(as.character(exp.table[,'mark'])), 'input')
    conditions <- unique(as.character(exp.table[,'condition']))
    if (length(conditions) < 2 & conf[['mode']] == 'differential') {
        stop("Mode 'differential' can only be used if two or more conditions are present.")
    }
    if (length(marks) < 2 & conf[['mode']] == 'combinatorial') {
        stop("Mode 'combinatorial' can only be used if two or more marks are present.")
    }
    
    ## Check if assembly must be present
    datafiles.clean <- sub('\\.gz$','', datafiles)
    format <- sapply(strsplit(datafiles.clean, '\\.'), function(x) { rev(x)[1] })
    if (any(format=='bed') & is.null(conf[['assembly']])) {
        stop("Please specify an 'assembly' for the BED files.")
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
    cat("- BROWSERFILES: Bed files for upload to the UCSC genome browser. It contains files with combinatorial states (*_combinations.bed.gz) and underlying peak calls (*_peaks.bed.gz). !!Always check the *_peaks.bed.gz files if you are satisfied with the peak calls. If not, there are ways to make the calls stricter (see section FAQ of the vignette).\n", file=savename, append=TRUE)
    cat("- -->combined<--: RData files with the combined results of the uni- and multivariate peak calling steps. This is what you want to use for downstream analyses. Contains combinedMultiHMM objects.\n", file=savename, append=TRUE)
    cat("    - combined_mode-separate.RData: Simple combination of peak calls (replicates considered) without multivariate analysis.\n", file=savename, append=TRUE)
    cat("    - combined_mode-combinatorial.RData: Combination of multivariate results for mode='combinatorial'.\n", file=savename, append=TRUE)
    cat("    - combined_mode-differential.RData: Combination of multivariate results for mode='differential'.\n", file=savename, append=TRUE)
    cat("    - combined_mode-full.RData: Combination of multivariate results for mode='full'.\n", file=savename, append=TRUE)
    cat("- multivariate: RData files with the results of the multivariate peak calling step. Contains multiHMM objects.\n", file=savename, append=TRUE)
    cat("- PLOTS: Several plots that are produced by default. Please check the plots in subfolder \"univariate-distributions\" for irregularities (see section \"Univariate Analysis\" of the vignette).\n", file=savename, append=TRUE)
    cat("- replicates: RData files with the result of the replicate peak calling step. Contains multiHMM objects.\n", file=savename, append=TRUE)
    cat("- univariate: RData files with the result of the univariate peak calling step. Contains uniHMM objects.\n", file=savename, append=TRUE)
    
    ## Write exclusive.table to file
    if (!is.null(conf[['exclusive.table']])) {
        if (is.character(conf[['exclusive.table']])) {
            conf[['exclusive.table']] <- utils::read.table(conf[['exclusive.table']], header=TRUE, comment.char='#')
        }
        check.exclusive.table(conf[['exclusive.table']])
        utils::write.table(conf[['exclusive.table']], file=file.path(outputfolder, 'exclusive_table.tsv'), col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
    }
    
    ## Make a copy of the conf file
    writeConfig(conf, configfile=file.path(outputfolder, 'chromstaR.config'))

    ## Write the experiment table to file
    utils::write.table(exp.table, file=file.path(outputfolder, 'experiment_table.tsv'), col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
    
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
  
  
    #==============
    ### Binning ###
    #==============
    messageU("Binning the data")
    ### Get chromosome lengths ###
    ## Get first bam file
    bamfile <- grep('bam$', datafiles, value=TRUE)[1]
    if (!is.na(bamfile)) {
        ptm <- startTimedMessage("Obtaining chromosome length information from file ", bamfile, " ...")
        chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
        stopTimedMessage(ptm)
    } else {
        ## Read chromosome length information
        if (is.character(conf[['assembly']])) {
            if (file.exists(conf[['assembly']])) {
                ptm <- startTimedMessage("Obtaining chromosome length information from file ", conf[['assembly']], " ...")
                df <- utils::read.table(conf[['assembly']], sep='\t', header=TRUE)
                stopTimedMessage(ptm)
            } else {
                ptm <- startTimedMessage("Obtaining chromosome length information from UCSC ...")
                df.chroms <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(conf[['assembly']])
                ## Get first bed file
                bedfile <- grep('bed$|bed.gz$', datafiles, value=TRUE)[1]
                if (!is.na(bedfile)) {
                    firstlines <- read.table(bedfile, nrows=10, skip=1) # skip 1 line in case of trackname
                    if (grepl('^chr',firstlines[1,1])) {
                        df <- df.chroms[,c('UCSC_seqlevel','UCSC_seqlength')]
                    } else {
                        df <- df.chroms[,c('NCBI_seqlevel','UCSC_seqlength')]
                    }
                }
                stopTimedMessage(ptm)
            }
        } else if (is.data.frame(conf[['assembly']])) {
            df <- conf[['assembly']]
        } else {
            stop("'assembly' must be either a data.frame with columns 'chromosome' and 'length' or a character specifying the assembly.")
        }
        chrom.lengths <- df[,2]
        names(chrom.lengths) <- df[,1]
        chrom.lengths <- chrom.lengths[!is.na(chrom.lengths) & !is.na(names(chrom.lengths))]
    }
    chrom.lengths.df <- data.frame(chromosome=names(chrom.lengths), length=chrom.lengths)
    
    ## Write chromosome length information to file
    utils::write.table(chrom.lengths.df, file=file.path(outputfolder, 'chrominfo.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    ### Make bins ###
    pre.bins.list <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=conf[['chromosomes']], binsizes=c(binsize, stepsize))
    pre.bins <- pre.bins.list[[1]]
    pre.bins.stepsize <- pre.bins.list[[2]]
    
    ### Count reads in bins ###
    if (!file.exists(binpath)) { dir.create(binpath) }
    parallel.helper <- function(ID, input, inputfile=NULL) {
        if (!input) {
            savename <- file.path(binpath, filenames[ID])
            savename.stepsize <- file.path(binpath, filenames.stepsize[ID])
            pairedEndReads <- exp.table[ID,'pairedEndReads']
            file <- datafiles[ID]
        } else {
            savename <- file.path(binpath, inputfilenames[basename(inputfile)])
            savename.stepsize <- file.path(binpath, inputfilenames.stepsize[basename(inputfile)])
            pairedEndReads <- exp.table[grep(basename(inputfile), exp.table$controlFiles),'pairedEndReads']
            if (any(pairedEndReads != pairedEndReads[1])) {
                stop("Multiple definitions of 'pairedEndReads' for file ", file, ".")
            }
            pairedEndReads <- pairedEndReads[1]
            file <- inputfile
        }
        if (!file.exists(savename)) {
            tC <- tryCatch({
                exp.table.input <- NULL
                if (!input) {
                    exp.table.input <- exp.table
                }
                binlist <- binReads(file=file, experiment.table=exp.table.input, ID=ID, assembly=chrom.lengths.df, pairedEndReads=pairedEndReads, binsizes=NULL, reads.per.bin=NULL, bins=list(pre.bins, pre.bins.stepsize), stepsizes=c(stepsize, stepsize), chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']])
                ptm <- startTimedMessage("Saving to file ", savename, " ...")
                bins <- binlist[[1]]
                save(bins, file=savename)
                stopTimedMessage(ptm)
                ptm <- startTimedMessage("Saving to file ", savename.stepsize, " ...")
                bins <- binlist[[2]]
                save(bins, file=savename.stepsize)
                stopTimedMessage(ptm)
            }, error = function(err) {
                stop(file,'\n',err)
            })
        }
    }
    ## Unparallelized
#     for (ID in IDs) {
#         parallel.helper(ID, input=FALSE)
#     }
#     for (file in inputfiles) {
#         parallel.helper(ID=NULL, input=TRUE, inputfile=file)
#     }
    ## Parallelized
    if (numcpu > 1) {
        ptm <- startTimedMessage("Binning data ...")
        temp <- foreach (ID = IDs, .packages=c("chromstaR")) %dopar% {
            parallel.helper(ID, input=FALSE)
        }
        stopTimedMessage(ptm)
    } else {
        for (ID in IDs) {
            parallel.helper(ID, input=FALSE)
        }
    }
    if (numcpu > 1) {
        ptm <- startTimedMessage("Binning input ...")
        temp <- foreach (file = inputfiles, .packages=c("chromstaR")) %dopar% {
            parallel.helper(ID=NULL, input=TRUE, inputfile=file)
        }
        stopTimedMessage(ptm)
    } else {
        for (file in inputfiles) {
            parallel.helper(ID=NULL, input=TRUE, inputfile=file)
        }
    }
    
    #==============================
    ### Univariate peak calling ###
    #==============================
    messageU("Calling univariate peaks")
    if (!file.exists(unipath)) { dir.create(unipath) }
    if (!file.exists(uniplotpath)) { dir.create(uniplotpath) }
    files <- file.path(binpath, filenames)
    names(files) <- names(filenames)
    files.stepsize <- file.path(binpath, filenames.stepsize)
    names(files.stepsize) <- names(filenames.stepsize)
    inputfiles.list <- lapply(strsplit(as.character(exp.table$controlFiles), '\\|'), function(x) { file.path(binpath, inputfilenames[x[!is.na(x)]]) })
    names(inputfiles.list) <- names(filenames)
    
    parallel.helper <- function(file, inputfiles) {
        ## Peak calling
        savename <- file.path(unipath, basename(file))
        if (!file.exists(savename)) {
            tC <- tryCatch({
                input.files <- NULL
                if (length(inputfiles)>0) {
                    input.files <- inputfiles
                }
                model <- callPeaksUnivariate(binned.data=file, input.data=input.files, eps=conf[['eps.univariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], read.cutoff.absolute=conf[['read.cutoff.absolute']], prefit.on.chr=conf[['prefit.on.chr']], keep.posteriors=FALSE, verbosity=0)
                # Replace counts.rpkm column with the binsize=stepsize counts, instead of the averaged ones
                binned.data.stepsize <- loadHmmsFromFiles(files = files.stepsize[names(file)], check.class = 'GRanges')[[1]]
                model$bins$counts.rpkm <- rpkm.vector(binned.data.stepsize$counts, binsize=stepsize)
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
        temp <- foreach (ID = names(files), .packages=c("chromstaR")) %dopar% {
            file <- files[ID]
            inputfiles <- inputfiles.list[[ID]]
            parallel.helper(file, inputfiles)
        }
        stopTimedMessage(ptm)
    } else {
        for (ID in names(files)) {
            file <- files[ID]
            inputfiles <- inputfiles.list[[ID]]
            parallel.helper(file, inputfiles)
        }
    }
  
  
    #================================
    ### Replicates                ###
    #================================
    if (mode == 'separate') {
      
        messageU("Analyzing replicates")
        if (!file.exists(reppath)) { dir.create(reppath) }
        
        files <- file.path(unipath, filenames)
        names(files) <- paste0(exp.table$mark, '-', exp.table$condition)
        markconditions <- unique(names(files))
        for (markcond in markconditions) {
            savename <- file.path(reppath, paste0(markcond, '_binsize', binsize.string, '_stepsize', stepsize.string, '.RData'))
            if (!file.exists(savename)) {
                mask <- names(files) == markcond
                repfiles <- files[mask]
                states <- stateBrewer(exp.table[mask,], mode='combinatorial')
                repmodel <- callPeaksMultivariate(repfiles, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']], temp.savedir=file.path(reppath, paste0(markcond, '_binsize', binsize.string, '_stepsize', stepsize.string, '_tempfiles')))
                ptm <- startTimedMessage("Saving to file ", savename, " ...")
                save(repmodel, file=savename)
                stopTimedMessage(ptm)
            } else {
                # repmodel <- loadHmmsFromFiles(savename, check.class=class.multivariate.hmm)[[1]]
            }
        }
        
        message("Combining replicate HMMs")
        modename <- 'separate'
        repfiles <- list.files(reppath, full.names=TRUE, pattern=paste0('binsize', binsize.string, '_stepsize', stepsize.string))
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
        tiles.per.cm <- 0.66
        min.tiles <- 5
        width <- max(min.tiles, length(combined.model$info$ID)) / tiles.per.cm + max(sapply(combined.model$info$ID, nchar)) / char.per.cm + legend.cm
        height <- max(min.tiles, length(combined.model$info$ID)) / tiles.per.cm + max(sapply(combined.model$info$ID, nchar)) / char.per.cm
        savename <- file.path(plotpath, 'read-count-correlation.pdf')
        ggplt <- heatmapCountCorrelation(combined.model, cluster=FALSE)
        ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
        savename <- file.path(plotpath, 'read-count-correlation-clustered.pdf')
        ggplt <- heatmapCountCorrelation(combined.model, cluster=TRUE)
        ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
        stopTimedMessage(ptm)
      
        #-------------------------
        ## Export browser files ##
        #-------------------------
        messageU("Exporting browser files")
        if (!file.exists(browserpath)) { dir.create(browserpath) }
        savename <- file.path(browserpath, paste0('combined_mode-', modename))
        trackname <- paste0('mode-', modename)
        if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
            exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='combinations', separate.files=TRUE)
        }
        if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
            exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='peaks', separate.files=TRUE)
        }
        invisible()
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
        char.per.cm <- 10
        legend.cm <- 3
        savename2 <- paste0(savename, '_transitionMatrix.pdf')
        if (!file.exists(savename2)) {
            ptm <- startTimedMessage("Making plots ...")
            multimodel <- suppressMessages( loadHmmsFromFiles(multimodel, check.class=class.multivariate.hmm)[[1]] )
            ggplt <- suppressMessages( heatmapTransitionProbs(multimodel) )
            width <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + legend.cm
            height <- length(levels(multimodel$bins$combination)) + max(sapply(levels(multimodel$bins$combination), nchar)) / char.per.cm + 1
            ggsave(savename2, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
            stopTimedMessage(ptm)
        }
    }
  
    ## Run multivariate depending on mode
    #--------------------
    if (mode == 'full') {
        multifile <- file.path(multipath, paste0('multivariate_mode-', mode, '.RData'))
        savename <- file.path(plotpath, paste0('multivariate_mode-', mode))
        if (!file.exists(multifile)) {
            files <- file.path(unipath, filenames)
            states <- stateBrewer(exp.table, mode=mode, exclusive.table=conf[['exclusive.table']])
            multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']], temp.savedir=file.path(multipath, paste0('multivariate_mode-', mode, '_tempfiles')))
            ptm <- startTimedMessage("Saving to file ", multifile, " ...")
            save(multimodel, file=multifile)
            stopTimedMessage(ptm)
            ## Plot transitions
            plothelper(savename, multimodel)
        } else {
            ## Plot transitions
            plothelper(savename, multifile)
        }
    
    #---------------------------
    } else if (mode == 'combinatorial') {
        for (condition in conditions) {
            messageU("condition = ", condition, underline='-', overline='-')
            multifile <- file.path(multipath, paste0('multivariate_mode-', mode, '_condition-', condition, '.RData'))
            savename <- file.path(plotpath, paste0('multivariate_mode-', mode, '_condition-', condition))
            if (!file.exists(multifile)) {
                mask <- exp.table[,'condition'] == condition
                files <- file.path(unipath, filenames)[mask]
                states <- stateBrewer(exp.table[mask,], mode=mode, exclusive.table=conf[['exclusive.table']])
                multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']], temp.savedir=file.path(multipath, paste0('multivariate_mode-', mode, '_condition-', condition, '_tempfiles')))
                ptm <- startTimedMessage("Saving to file ", multifile, " ...")
                save(multimodel, file=multifile)
                stopTimedMessage(ptm)
                ## Plot transitions
                plothelper(savename, multimodel)
            } else {
                ## Plot transitions
                plothelper(savename, multifile)
            }
        }
    
    #--------------------------------
    } else if (mode == 'differential') {
        for (mark in marks) {
            messageU("mark = ", mark, underline='-', overline='-')
            multifile <- file.path(multipath, paste0('multivariate_mode-', mode, '_mark-', mark, '.RData'))
            savename <- file.path(plotpath, paste0('multivariate_mode-', mode, '_mark-', mark))
            if (!file.exists(multifile)) {
                mask <- exp.table[,'mark'] == mark
                files <- file.path(unipath, filenames)[mask]
                states <- stateBrewer(exp.table[mask,], mode=mode, exclusive.table=conf[['exclusive.table']])
                multimodel <- callPeaksMultivariate(files, use.states=states, max.states=conf[['max.states']], eps=conf[['eps.multivariate']], max.iter=conf[['max.iter']], max.time=conf[['max.time']], num.threads=conf[['numCPU']], per.chrom=conf[['per.chrom']], keep.posteriors=conf[['keep.posteriors']], temp.savedir=file.path(multipath, paste0('multivariate_mode-', mode, '_mark-', mark, '_tempfiles')))
                ptm <- startTimedMessage("Saving to file ", multifile, " ...")
                save(multimodel, file=multifile)
                stopTimedMessage(ptm)
                ## Plot transitions
                plothelper(savename, multimodel)
            } else {
                ## Plot transitions
                plothelper(savename, multifile)
            }
        }
    }

  
    #================================
    ## Combine multiple conditions ##
    #================================
    messageU("Combining multivariate HMMs")
    multifiles <- list.files(multipath, full.names=TRUE, pattern=paste0('mode-',mode))

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
    ## Plot correlations
    ptm <- startTimedMessage("Plotting read count correlation ...")
    char.per.cm <- 10
    legend.cm <- 3
    tiles.per.cm <- 0.66
    min.tiles <- 5
    width <- max(min.tiles, length(combined.model$info$ID)) / tiles.per.cm + max(sapply(combined.model$info$ID, nchar)) / char.per.cm + legend.cm
    height <- max(min.tiles, length(combined.model$info$ID)) / tiles.per.cm + max(sapply(combined.model$info$ID, nchar)) / char.per.cm
    savename <- file.path(plotpath, 'read-count-correlation.pdf')
    ggplt <- heatmapCountCorrelation(combined.model, cluster=FALSE)
    ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
    savename <- file.path(plotpath, 'read-count-correlation-clustered.pdf')
    ggplt <- heatmapCountCorrelation(combined.model, cluster=TRUE)
    ggsave(savename, plot=ggplt, width=width, height=height, limitsize=FALSE, units='cm')
    stopTimedMessage(ptm)
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    messageU("Exporting browser files")
    if (!file.exists(browserpath)) { dir.create(browserpath) }
    savename <- file.path(browserpath, paste0('combined_mode-', mode))
    trackname <- paste0('mode-', mode)
    if (!file.exists(paste0(savename, '_combinations.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='combinations', separate.files=TRUE)
    }
    if (!file.exists(paste0(savename, '_peaks.bed.gz'))) {
        exportCombinedMultivariate(combined.model, filename=savename, trackname=trackname, what='peaks', separate.files=TRUE)
    }
  
    total.time <- proc.time() - total.time
    message("==> Total time spent: ", round(total.time[3]), "s <==")
  
  
}
