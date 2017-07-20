

#' Read chromstaR configuration file
#'
#' Read a chromstaR configuration file into a list structure. The configuration file has to be specified in INI format. R expressions can be used and will be evaluated.
#'
#' @param configfile Path to the configuration file
#' @return A \code{list} with one entry for each element in \code{configfile}.
#' @author Aaron Taudt
#' @importFrom utils read.table
readConfig <- function(configfile) {

    connection <- file(configfile) 
    Lines  <- readLines(connection) 
    close(connection) 

    Lines <- chartr("[]", "==", Lines) # change section headers 
    Lines <- gsub(" ", "", Lines) # no spaces

    connection <- textConnection(Lines) 
    data <- utils::read.table(connection, as.is = TRUE, sep = "=", fill = TRUE, quote="") 
    close(connection) 
    names(data) <- c('argument','value','section')

    L <- data$argument == "" # location of section breaks 
    data$section <- data$value[which(L)[cumsum(L)]]
    data <- data[data$argument!="",]

    configlist <- list() 
    ToParse <- paste0("configlist <- list(", paste(paste(data$argument, data$value, sep="="), collapse=', '), ")")

    eval(parse(text=ToParse)) 

    return(configlist) 
} 

#' Write chromstaR configuration file
#'
#' Write a chromstaR configuration file from a list structure.
#'
#' @param conf A list structure with parameter values. Each entry will be written in one line.
#' @param configfile Filename of the outputfile.
#' @return \code{NULL}
#' @author Aaron Taudt
#' @importFrom utils write.table
writeConfig <- function(conf, configfile) {

    ## Printing function
    formatstring <- function(string) {
        if (is.factor(string)) {
            string <- as.character(string)
        }
        if (is.character(string) & length(string)>1) {
            string <- paste0("c('",paste0(string,collapse="','"),"')")
        } else if (is.character(string) & length(string)==1) {
            string <- paste0("'",string,"'")
        } else if (is.numeric(string) & length(string)>1) {
            string <- paste0("c(",paste0(string,collapse=','),")")
        } else if (is.numeric(string) & length(string)==1) {
            string <- string
        } else if (is.null(string)) {
            string <- "NULL"
        } else if (is.data.frame(string)) {
            if (all(names(string) %in% c('chromosome','length'))) {
                string <- paste0("'", file.path(dirname(configfile), 'chrominfo.tsv'), "'")
            } else if (all(names(string) %in% c('mark','group'))) {
                string <- paste0("'", file.path(dirname(configfile), 'exclusive_table.tsv'), "'")
            }
        }
        return(string)
    }
        
    f <- file(configfile, open='w')
    cat("#============== chromstaR configuration file ===============#\n", file=f)
    cat("\n[General]\n", file=f)
    for (i1 in c('numCPU')) {
        cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
    }
    cat("\n[Binning]\n", file=f)
    for (i1 in c('binsize', 'stepsize', 'assembly', 'chromosomes', 'remove.duplicate.reads', 'min.mapq')) {
        cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
    }
    cat("\n[Univariate]\n", file=f)
    for (i1 in c('prefit.on.chr', 'eps.univariate', 'max.time', 'max.iter', 'read.cutoff.absolute')) {
        cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
    }
    cat("\n[Multivariate]\n", file=f)
    for (i1 in c('mode', 'eps.multivariate', 'max.states', 'per.chrom', 'keep.posteriors', 'exclusive.table')) {
        cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
    }
    close(f, type='w')
}
