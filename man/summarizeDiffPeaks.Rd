\name{summarizeDiffPeaks}
\alias{summarizeDiffPeaks}
\title{
Summary information about differential peak calls
}
\description{
Get information about the differential peak calls in a \code{\link{chromstaR_multivariateHMM}} object.
}
\usage{
summarizeDiffPeaks(multi.hmm)
}
\arguments{
\item{multi.hmm}{A \code{\link{chromstaR_multivariateHMM}} object.}
}

\value{
A list object with various entries.
}
\keyword{summary}
\author{Aaron Taudt}
\examples{
## Load example multivariate Hidden Markov Model
data(example.multi.HMM)
## Summary information about differentially modified peaks
summarizeDiffPeaks(example.multi.HMM)
}
