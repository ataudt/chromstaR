\name{Convert to GenomicRanges}
\alias{bed2GRanges}

\title{
Convert contents of a BED file to \link{GenomicRanges} objects.
}
\description{
Convert contents of a BED file to GRanges objects. Only coordinates and information per position (e.g. reads, posteriors) are converted, other information specific to the HMM (e.g. loglikelihood, model parameters) are not converted.
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
}
