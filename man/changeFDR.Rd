\name{changeFDR}
\alias{changeFDR}
\title{
Change the false discovery rate of a Hidden Markov Model
}
\description{
Adjusts the peak calls of a \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} object with the given false discovery rate (FDR).
}
\usage{
changeFDR(model, FDR=0.5, separate.zeroinflation=TRUE)
}
\arguments{
\item{model}{A \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} object with posteriors.}
\item{FDR}{False discovery rate.}
\item{separate.zeroinflation}{Only for \code{\link{chromstaR_univariateHMM}} objects: If set to TRUE, state 'zero-inflation' will be treated separately, otherwise it will be merged with state 'unmodified'.}
}

\details{
Peaks are called if the posterior for a state (univariate) or sample (multivariate) are \code{>= 1-FDR}.
}
\value{
The input object is returned with adjusted peak calls.
}
\keyword{FDR}
\keyword{fdr}
\author{Aaron Taudt}
\examples{
}
