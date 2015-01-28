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
}

