\name{binned.data.example}
\docType{data}
\alias{binned.data.example.1}
\alias{binned.data.example.2}

\title{Example datasets for demonstrating the workflow of a \code{\link{chromStar}} script}
\description{These datasets contain aligned reads for chromosome 1 and two histone modifications (H3k27me3 and H3k36me3), binned into bins of size 2000.}
\usage{
data(binned.data.example.1)
data(binned.data.example.2)
}
\format{
A data.frame with 4 columns. The first 3 columns contain genomic coordinates ('chromosome', 'start' and 'end') defining a bin. The 4th column contains the aligned reads in that bin.
}
\source{
Wiggler: Normalized genome-wide signal coverage tracks of ENCODE data
\url{http://www.broadinstitute.org/~anshul/projects/encode/rawdata/mapped/jan2011/noMultiMapTagAlign/wgEncodeBroadHistone/female/}
}

\keyword{binned}
\keyword{data}
