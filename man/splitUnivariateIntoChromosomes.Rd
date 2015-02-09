\name{splitUnivariateIntoChromosomes}
\alias{splitUnivariateIntoChromosomes}

\title{
Split a \code{\link{chromstaR_univariateHMM}} into chromosomes
}
\description{
Split a \code{\link{chromstaR_univariateHMM}} object into chromosomes to enable parallelized fitting of the multivariate Hidden Markov Model on separate chromosomes instead of fitting on the full genome. This is mainly done for performance reasons.
}
\usage{
splitUnivariateIntoChromosomes(uni.hmm, filename=NULL)
}
\arguments{
\item{uni.hmm}{A \code{\link{chromstaR_univariateHMM}} object.}
\item{filename}{The file name where the split objects will be stored. Each chromosome is saved in a separate file, with the chromosome name appended to the file name. If \code{filename} is not specified, a list of \code{\link{chromstaR_univariateHMM}}s is returned.}
}

\value{
A list of \code{\link{chromstaR_univariateHMM}} objects or NULL, depending on option \code{filename}.
}
\keyword{split}
\author{Aaron Taudt}
\seealso{\code{\link{mergeChromsFromMultivariates}}}
\examples{
## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue (chr22)
bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz", package="chromstaR")
## Bin the BED file into bin size 800bp
binned.data <- bed2binned(bedfile, assembly='hg19', binsize=800, save.as.RData=FALSE)
## Fit the univariate Hidden Markov Model (HMM)
hmm <- callPeaksUnivariate(binned.data, ID='example_H3K36me3', max.time=60)
## Check if the fit is ok
plot(hmm, type='histogram')
## Split the fitted HMM into separate chromosomes (not sensible here because the example dataset consists only of 1 chromosome)
\donttest{splitUnivariateIntoChromosomes(hmm, filename='chromstaR-example_H3K36me3')}
}

