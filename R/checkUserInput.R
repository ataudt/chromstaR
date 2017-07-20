check.positive.integer <- function(testvar) {
    if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
    if (length(testvar)>1) return(2)
    if (testvar <= 0) return(3)
    if (testvar != as.integer(testvar)) return(4)
    return(0)
}
check.nonnegative.integer <- function(testvar) {
    if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
    if (length(testvar)>1) return(2)
    if (testvar < 0) return(3)
    if (testvar != as.integer(testvar)) return(4)
    return(0)
}
check.positive.integer.vector <- function(testvec) {
    if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
    for (elem in testvec) {
        if (elem <= 0) return(2)
        if (elem != as.integer(elem)) return(3)
    }
    return(0)
}
check.nonnegative.integer.vector <- function(testvec) {
    if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
    for (elem in testvec) {
        if (elem < 0) return(2)
        if (elem != as.integer(elem)) return(3)
    }
    return(0)
}
check.positive <- function(testvar) {
    if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
    if (length(testvar)>1) return(2)
    if (testvar <= 0) return(3)
    return(0)
}
check.positive.vector <- function(testvec) {
    if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
    for (elem in testvec) {
        if (elem <= 0) return(2)
    }
    return(0)
}
check.nonnegative.vector <- function(testvec) {
    if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
    for (elem in testvec) {
        if (elem < 0) return(2)
    }
    return(0)
}
check.integer <- function(testvar) {
    if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
    if (length(testvar)>1) return(2)
    if (testvar != as.integer(testvar)) return(3)
    return(0)
}

check.univariate.modellist <- function(modellist) {
    if (!is(modellist,"list")) return(1)
    for (model in modellist) {
        if (!is(model,class.univariate.hmm)) return(2)
    }
    return(0)
}
check.univariate.model <- function(model) {
    if (!is(model,class.univariate.hmm)) return(1)
    return(0)
}

check.multivariate.modellist <- function(modellist) {
    if (!is(modellist,"list")) return(1)
    for (model in modellist) {
        if (!is(model,class.multivariate.hmm)) return(2)
    }
    return(0)
}
check.multivariate.model <- function(model) {
    if (!is(model,class.multivariate.hmm)) return(1)
    return(0)
}

check.logical <- function(testbool) {
    if (!is(testbool,"logical")) return(1)
    if (length(testbool)>1) return(2)
    return(0)
}

check.experiment.table <- function(experiment.table) {
    err <- 0
    if (!is.data.frame(experiment.table)) err <- 1
    if (any(names(experiment.table) != c('file', 'mark', 'condition', 'replicate', 'pairedEndReads', 'controlFiles'))) err <- 2
    if (err == 1 | err == 2) {
        stop("Argument 'experiment.table' expects a data.frame with columns 'file', 'mark', 'condition', 'replicate', 'pairedEndReads' and 'controlFiles'.")
    }
    conditions <- unique(experiment.table$condition)
    if (any(grepl('[[:punct:]]', conditions)) | any(grepl('[[:punct:]]', conditions)) | any(grepl('^[0-9]', conditions))) {
        stop("Column 'condition' of the experiment.table cannot contain special characters or spaces or start with a number.")
    }
    marks <- unique(experiment.table$mark)
    if (any(grepl('[[:punct:]]', marks)) | any(grepl('[[:punct:]]', marks))) {
        stop("Column 'mark' of the experiment.table cannot contain special characters or spaces.")
    }
    replicates <- unique(experiment.table$replicate)
    if (any(grepl('[[:punct:]]', replicates)) | any(grepl('[[:punct:]]', replicates))) {
        stop("Column 'replicate' of the experiment.table cannot contain special characters or spaces.")
    }
    IDs <- paste0(experiment.table$mark, '-', experiment.table$condition, '-rep', experiment.table$replicate)
    tab <- table(IDs)
    if (any(duplicated(IDs))) {
        stop("Duplicated IDs detected. Check your experiment.table: ", paste(paste0(names(tab), ' (', tab, ')'), collapse = ', '))
    }
    
    return(err)
}

check.exclusive.table <- function(exclusive.table) {
    err <- 0
    if (!is.data.frame(exclusive.table)) err <- 1
    if (any(names(exclusive.table) != c('mark', 'group'))) err <- 2
    if (err > 0) {
        stop("Argument 'exclusive.table' expects a data.frame with columns 'mark', 'group'.")
    }
    return(err)
}

