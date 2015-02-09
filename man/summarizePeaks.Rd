\name{summarizePeaks}
\alias{summarizePeaks}
\title{
Summary information about called peaks
}
\description{
Get information about the peak calls in a \code{\link{chromstaR_multivariateHMM}} object.
}
\usage{
summarizePeaks(multi.hmm)
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
## Summary information about peaks
summarizePeaks(example.multi.HMM)
}
