\name{binned.data}
\docType{data}
\alias{binned.data}

\title{
Format of the data.frame containing binned read counts
}
\description{
This data.frame is ouput of the binning functions \code{\link{bam2binned}}, \code{\link{bed2binned}} and \code{\link{bedGraph2binned}}. The read counts can be visualized in a Genome Browser using the export function \code{\link{binned2wiggle}}. This data.frame is input for the univariate region calling with function \code{\link{univariate.from.binned.data}}.
}
\format{
A data.frame with 4 columns. Columns 1 to 3 are chromosomal coordinates and column 4 is the read count in the corresponding bin. Example:
\tabular{rrrr}{
chrom \tab start \tab end \tab reads \cr
... \tab ... \tab ... \tab ... \cr
chr1 \tab 28000 \tab 28199 \tab 0 \cr
chr1 \tab 28200 \tab 28399 \tab 0 \cr
chr1 \tab 28400 \tab 28599 \tab 3 \cr
chr1 \tab 28600 \tab 28799 \tab 0 \cr
... \tab ... \tab ... \tab ... \cr
}
}

\keyword{binned}
\keyword{data}
\keyword{binned.data}
\seealso{\code{\link{binned2GRanges}}, \code{\link{binned2wiggle}}}
