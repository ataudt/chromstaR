\name{bed2GRanges}
\alias{bed2GRanges}

\title{
Read a BED file into a \link{GRanges} objects.
}
\description{
Read a BED file into a \link{GRanges} object. Only coordinates are read.
}
\usage{
bed2GRanges(bedfile, chrom.length.file, skip=1, binsize=NULL)
}
\arguments{
\item{bedfile}{A BED file.}
\item{chrom.length.file}{
A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.}
\item{skip}{Number of lines to skip in BED file.}
\item{binsize}{Force the resulting \link{GRanges} to have this bin size.}
}

\value{
A \link{GRanges} object.
}
\keyword{GRanges}
\author{Aaron Taudt}
\examples{
## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue (chr22)
bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.chr22.bed.gz", package="chromstaR")
## Read as GRanges
gr <- bed2GRanges(bedfile)
}
