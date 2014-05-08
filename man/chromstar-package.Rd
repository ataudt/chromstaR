\name{chromstar}
\alias{chromstar-package}
\alias{chromstar}
\docType{package}
\title{
Combinatorial analysis of ChIP-seq histone modification data using Hidden Markov Models
}
\description{
\strong{chromstar} is an R package for analysis of histone modification data (ChIP-seq). Both univariate and multivariate analysis of multiple histone modifications are possible. A Hidden Markov Model with a mixture of Negative Binomials as emission densities is used to classify the data. Please read our publication for a detailed description of the method.
TODO: insert publication
}
\details{
\tabular{ll}{
Package: \tab chromstar\cr
Type: \tab Package\cr
Version: \tab 0.8\cr
Date: \tab 2014-01-23\cr
Depends: \tab R (>= 3.0.1)\cr
Imports: \tab gsl\cr
Suggests: \tab ggplot2\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Packaged: \tab 2013-06-20 14:51:55 UTC; heinig\cr
Built: \tab R 3.0.2; x86_64-pc-linux-gnu; 2014-01-31 10:10:41 UTC; unix\cr
}
}
\section{Index}{
\preformatted{
binned.data.example.1   Example datasets for demonstrating the workflow
                        of a \strong{chromstar} script
collapse.bins           Collapse consecutive bins
combinatorial.states    Get the (decimal) combinatorial states of a
                        list of HMM models
histoneHMM-package      Analysis of histone modification data using a
                        HMM
make.ucsc.file.combinatorial
                        Convert univariate HMM models for upload to the
                        UCSC genome browser
make.ucsc.file.multivariate
                        Convert a multivariate HMM model for upload to
                        the UCSC genome browser
make.ucsc.file.univariate
                        Convert univariate HMM models for upload to the
                        UCSC genome browser
multivariate.from.univariate.results
                        Perform the multivariate HMM with univariate
                        results
plot.distribution       Plot the emission distributions that were
                        fitted to the read data
univariate.from.binned.data
                        Perform the univariate HMM with binned data
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Maria Colome Tatche, Aaron Taudt, Matthias Heinig 

Maintainer: Aaron Taudt <a.s.taudt@umcg.nl>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
