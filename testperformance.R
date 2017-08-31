library(devtools)
library(chromstaR)

load_all()

## Load simulated ground truth ##
load('~/Desktop/Arbeit/DIR/chromstaR/simulated_counts.RData')
states <- dec2bin(sim$bins$state)
counts <- sim$bins$counts
exp <- data.frame(file=file.path('data', paste0(names(sim$reads), '.bed.gz')), mark='H3K4me3', condition=paste0('SIM', 1:4), replicate=1, pairedEndReads=FALSE, controlFiles=NA)
exp$ID <- paste0(exp$mark, '-', exp$condition, '-rep', exp$replicate)


## Make GRanges for usage with chromstaR ##
bins <- GRanges(seqnames='chr12', ranges = IRanges(start=1:nrow(counts), end=2:(nrow(counts)+1)))

## Univariate fits ##
models <- list()
for (i1 in 1:ncol(counts)) {
  binned.data <- bins
  binned.data$counts <- matrix(counts[,i1], nrow=nrow(counts))
  colnames(binned.data$counts) <- '0'
  attr(binned.data, 'info') <- exp[i1,]
  models[[i1]] <- callPeaksUnivariate(binned.data, max.time=600, eps=1e-4, max.iter=1)
  rownames(sim$emissions[[i1]]) <- rownames(models[[i1]]$distributions)
  models[[i1]]$distributions <- sim$emissions[[i1]]
}

model.copula <- callPeaksMultivariate(models, use.states = NULL, eps=1e-2, max.iter=3)
model.influence <- callPeaksInfluence(models, eps=1e-2, max.iter=3)

states.copula <- dec2bin(model.copula$bins$state, colnames = model.copula$info$ID)
states.influence <- dec2bin(model.influence$bins$state, colnames = model.influence$info$ID)

## Table for performance metrics
df <- data.frame(truth=as.vector(states), prediction=as.vector(states.copula))
table(df)

## Differential states
diff.states <- rowSums(states) > 0 & rowSums(states) < 4

