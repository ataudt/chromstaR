\name{chromstaR}
\alias{chromstaR-package}
\alias{chromstaR}
\docType{package}
\title{
Combinatorial analysis of ChIP-seq histone modification data using Hidden Markov Models
}
\description{
\pkg{chromstaR} is an R package for the combinatorial or differential analysis of multiple ChIP-seq experiments. It was developed for histone modifications with a broad profile but is also suitable for the analysis of transcription factor binding data. A Hidden Markov Model with a mixture of Negative Binomials as emission densities is used to call peaks. Please read our publication for a detailed description of the method.
TODO: insert publication
}
\details{
\tabular{ll}{
Package: \tab chromstaR\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2015-01-27\cr
Depends: \tab R (>= 3.1.0), GenomicRanges, IRanges\cr
Imports: \tab Rsamtools, BiocGenerics, GenomicAlignments\cr
Suggests: \tab ggplot2, reshape2, grid, rtracklayer, mvtnorm, biomaRt\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Packaged: \tab 2013-06-20 14:51:55 UTC; Aaron Taudt\cr
Built: \tab R 3.0.2; x86_64-pc-linux-gnu; 2014-01-31 10:10:41 UTC; unix\cr
}
}
\section{Index}{
\preformatted{
binned.data.example.1   Example datasets for demonstrating the workflow
                        of a \pkg{chromstaR} script
collapseBins           Collapse consecutive bins
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
callPeaksUnivariate
                        Perform the univariate HMM with binned data
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Aaron Taudt, Maria Colome Tatche, Matthias Heinig 

Maintainer: Aaron Taudt <a.s.taudt@umcg.nl>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
\code{\link{GRanges}}
}
\examples{
}
