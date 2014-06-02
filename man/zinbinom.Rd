\name{Zero-inflated NegBinomial}
\alias{dzinbinom}
\alias{pzinbinom}
\alias{qzinbinom}
\alias{rzinbinom}

\title{
The zero-inflated Negative Binomial Distribution
}
\description{
Density, distribution function, quantile function and random generation for the zero-inflated negative binomial distribution with parameters \option{w}, \option{size} and \option{prob}.
}
\usage{
dzinbinom(x, w, size, prob, mu)
pzinbinom(q, w, size, prob, mu, lower.tail=TRUE)
qzinbinom(p, w, size, prob, mu, lower.tail=TRUE)
rzinbinom(n, w, size, prob, mu)
}
\arguments{
Please refer to \code{\link{dnbinom}} for the full list of arguments. The only additional parameter for the zero-inflated Negative Binomial is
\item{w}{Weight of the zero-inflation (between 0 and 1).}
}
\details{
TODO: write down formula
}
\value{
\option{dnbinom} gives the density, \option{pnbinom} gives the distribution function, \option{qnbinom} gives the quantile function, and \option{rnbinom} generates random deviates.
}
\source{}
\seealso{\code{\link{dnbinom}}}
\keyword{zero-inflation}
\keyword{inflation}
\keyword{inflated}
