\name{unis2pseudomulti}
\alias{unis2pseudomulti}

\title{
Combine univariate HMMs to a multivariate HMM
}
\description{
Combine multiple \code{\link{chromstaR_univariateHMM}}s to a \code{\link{chromstaR_multivariateHMM}} without running \code{\link{callPeaksMultivariate}}. This should only be done for comparison purposes.
}
\usage{
unis2pseudomulti(uni.hmm.list)
}
\arguments{
\item{uni.hmm.list}{A list of \code{\link{chromstaR_univariateHMM}}s, e.g. \code{list(hmm1, hmm2, ...)}.}
}

\details{
Use this function if you want to combine ChIP-seq samples without actually running a multivariate Hidden Markov Model. The resulting object will be of class \code{\link{chromstaR_multivariateHMM}} but will not be truly multivariate.
}
\value{
A \code{\link{chromstaR_multivariateHMM}} object.
}
\author{Aaron Taudt}
\examples{
## Get example BED-files with ChIP-seq reads for H3K36me3 in 7 different brain tissues (chr22)
bedfiles <- list.files(system.file("extdata/brain/", package="chromstaR"), full=T)
## Bin the data into bin size 800bp and build the univariate Hidden Markov Model (HMM)
uni.HMMs <- list()
for (bedfile in bedfiles) {
  binned.data <- bed2binned(bedfile, assembly='hg19', binsize=800, save.as.RData=F)
  uni.HMMs[[bedfile]] <- callPeaksUnivariate(binned.data, ID=basename(bedfile), max.time=30, eps=0.01)
}
## Combine the univariate HMMs without fitting a multivariate HMM
pseudo.multi.HMM <- unis2pseudomulti(uni.HMMs)
}

