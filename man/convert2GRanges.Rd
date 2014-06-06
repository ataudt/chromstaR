\name{Convert to GenomicRanges}
\alias{hmm2GRanges}
\alias{binned2GRanges}

\title{
Convert various \code{\link{chromstar}} objects into \code{\link{GenomicRanges}} objects.
}
\description{
Convert \code{\link{binned.data}}, \code{\link{uni.hmm}}s and \code{\link{multi.hmm}}s to GRanges objects. Only coordinates and information per position (e.g. reads, posteriors) are converted, other information specific to the HMM (e.g. loglikelihood, model parameters) are not converted.
}
\usage{
binned2GRanges(binned.data)

hmm2GRanges(hmm, reduce=TRUE)
}
\arguments{
\item{binned.data}{A data.frame with coordinates and binned read counts (see \code{\link{binned.data}}).}
\item{hmm}{Either a \code{\link{uni.hmm}} or \code{\link{multi.hmm}} object, produced by \code{\link{univariate.from.binned.data}} or \code{\link{multivariate.from.univariate.results}}, respectively.}
\item{reduce}{Either TRUE or FALSE, indicating whether consecutive bins with the same state should be merged into bigger regions.}
}

\value{
A GRanges object is returned.
}
\keyword{GRanges}
\author{Aaron Taudt}
\examples{
}
