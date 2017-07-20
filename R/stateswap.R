stateswap <- function(model) {
  
    distr <- model$distributions
    mean1 <- distr['unmodified','mu']
    mean2 <- distr['modified','mu']
    valat0_1 <- dnbinom(x = 0, size = distr['unmodified','size'], prob = distr['unmodified','prob'])
    valat0_2 <- dnbinom(x = 0, size = distr['modified','size'], prob = distr['modified','prob'])
    maxcount <- max(model$bincounts$counts[,'0'])
    valatmax_1 <- dnbinom(x = maxcount, size = distr['unmodified','size'], prob = distr['unmodified','prob'])
    valatmax_2 <- dnbinom(x = maxcount, size = distr['modified','size'], prob = distr['modified','prob'])
    
    meanL <- mean2 >= mean1
    valat0L <- valat0_1 >= valat0_2
    valatmaxL <- valatmax_2 >= valatmax_1
    
    if (sum(c(mean=meanL, valat0=valat0L, valatmax=valatmaxL)) <= 1) {
        ## Swap states ##
        model$distributions <- model$distributions[c(1,3,2),]
        rownames(model$distributions) <- rownames(model$distributions)[c(1,3,2)]
        state <- model$bins$state
        state[model$bins$state == 'unmodified'] <- 'modified'
        state[model$bins$state == 'modified'] <- 'unmodified'
        model$bins$state <- state
        model$weights <- model$weights[c(1,3,2)]
        names(model$weights) <- names(model$weights)[c(1,3,2)]
        trans <- model$transitionProbs
        trans[1,2] <- model$transitionProbs[1,3]
        trans[1,3] <- model$transitionProbs[1,2]
        trans[2,1] <- model$transitionProbs[3,1]
        trans[3,1] <- model$transitionProbs[2,1]
        trans[2,2] <- model$transitionProbs[3,3]
        trans[3,3] <- model$transitionProbs[2,2]
        trans[2,3] <- model$transitionProbs[3,2]
        trans[3,2] <- model$transitionProbs[2,3]
        model$transitionProbs <- trans
        model$startProbs <- model$startProbs[c(1,3,2)]
        names(model$startProbs) <- names(model$startProbs)[c(1,3,2)]
        model$bins$posteriors <- model$bins$posteriors[,c(1,3,2)]
        colnames(model$bins$posteriors) <- colnames(model$bins$posteriors)[c(1,3,2)]
        model$bins$posterior.modified <- model$bins$posteriors[,'modified']
    }
    
    ## Check if swap was useful and issue warning if not ##
    distr <- model$distributions
    mean1 <- distr['unmodified','mu']
    mean2 <- distr['modified','mu']
    valat0_1 <- dnbinom(x = 0, size = distr['unmodified','size'], prob = distr['unmodified','prob'])
    valat0_2 <- dnbinom(x = 0, size = distr['modified','size'], prob = distr['modified','prob'])
    maxcount <- max(model$bincounts$counts[,'0'])
    valatmax_1 <- dnbinom(x = maxcount, size = distr['unmodified','size'], prob = distr['unmodified','prob'])
    valatmax_2 <- dnbinom(x = maxcount, size = distr['modified','size'], prob = distr['modified','prob'])
    
    meanL <- mean2 >= mean1
    valat0L <- valat0_1 >= valat0_2
    valatmaxL <- valatmax_2 >= valatmax_1
    
    if (sum(c(mean=meanL, valat0=valat0L, valatmax=valatmaxL)) <= 1) {
        warning(model$info$ID, ": State assignment is messed up. Please check the fits in 'outputfolder/PLOTS/univariate-distributions'. Possible fixes are 1) pick another binsize or 2) downsample your data or 3) use data of better quality.")
    }
    
    
    return(model)
}
