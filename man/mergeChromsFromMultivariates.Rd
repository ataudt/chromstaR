\name{mergeChromsFromMultivariates}
\alias{mergeChromsFromMultivariates}

\title{
Merge several \code{\link{chromstaR_multivariateHMM}}s into one object
}
\description{
Merge several \code{\link{chromstaR_multivariateHMM}}s into one object. This can be done to merge fits for separate chromosomes into one object for easier handling.
}
\usage{
mergeChromsFromMultivariates(multi.hmm.list, filename=NULL)
}
\arguments{
\item{multi.hmm.list}{A list of \code{\link{chromstaR_multivariateHMM}} objects.}
\item{filename}{The file name where the merged object will be stored. If \code{filename} is not specified, a \code{\link{chromstaR_multivariateHMM}} is returned.}
}

\value{
A \code{\link{chromstaR_multivariateHMM}} object or NULL, depending on option \code{filename}.
}
\keyword{merge}
\author{Aaron Taudt}
\seealso{\code{\link{splitUnivariateIntoChromosomes}}}
\examples{
}

