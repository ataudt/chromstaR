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
}

