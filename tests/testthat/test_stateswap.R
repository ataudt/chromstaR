messageU("State swap")

### Test all binned data in folder stateswap for correct convergence ###
files <- list.files('stateswap', full.names=TRUE)

## H3K27Ac_GMP_Rep1
ID <- "H3K27Ac_GMP_Rep1"
file <- grep(ID, files, value=T)
unimodel <- callPeaksUnivariate(file, eps=1)
expect_that(unimodel$weights['modified'], is_less_than(0.2))
