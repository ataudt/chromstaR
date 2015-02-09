\name{plotting functions for chromstaR}
\alias{plot}
\alias{plot.chromstaR_univariateHMM}
\alias{plot.chromstaR_multivariateHMM}

\title{
Plotting functions for the \pkg{\link{chromstaR}} package
}
\description{
Plots for \code{\link{chromstaR_univariateHMM}} and \code{\link{chromstaR_multivariateHMM}} objects.
}
\usage{
\method{plot}{chromstaR_univariateHMM}(x, type, \dots)
\method{plot}{chromstaR_multivariateHMM}(x, type, \dots)
}
\arguments{
\item{x}{
A \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} object.
}
\item{type}{Type of the plot. The following types are accepted for \code{\link{chromstaR_univariateHMM}} plots:
	\describe{
		\item{\code{histogram}}{A histogram with fitted mixture distribution.}
		\item{\code{boxplot}}{A boxplot of read counts for the different states.}
		\item{\code{normalTransformation}}{A histogram of transformed read counts.}
	}
	The following types are accepted for \code{\link{chromstaR_multivariateHMM}} objects:
	\describe{
		\item{\code{transitionMatrix}}{A heatmap with entries of the transition matrix.}
		\item{\code{histograms}}{Fitted histograms of all underlying univariate distributions.}
	}
}
\item{\dots}{Additional arguments that can be passed to the different plot-types:
	\describe{
		\item{\code{state}}{The modification state \code{c('unmodified','modified')} for \code{type=c('histogram','normalTransformation')}.}
		\item{\code{chrom}}{Chromosome to plot for \code{type=c('histogram')}.}
		\item{\code{start}}{Start position to plot for \code{type=c('histogram')}.}
		\item{\code{end}}{End position to plot for \code{type=c('histogram')}.}
	}
}
}

\value{
A \link[ggplot2:ggplot]{ggplot} object or a list of such objects.
}
\keyword{histogram}
\keyword{distributions}
\author{Aaron Taudt}

