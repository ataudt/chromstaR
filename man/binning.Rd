\name{Convert to binned data}
\alias{bam2binned}
\alias{bed2binned}
\alias{bedGraph2binned}

\title{
Convert aligned reads from various file formats into read counts in equidistant windows (bins)
}
\description{
Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
}
\usage{
bam2binned(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE)

bed2binned(bedfile, chrom.length.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE)

bedGraph2binned(bedGraphfile, chrom.length.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE)
}
\arguments{
\item{bamfile}{A file with aligned reads in .bam format.}
\item{bamindex}{A bam-index file (.bai). If you miss this file, generate it with the command 'samtools index bamfile.bam' of the SAMtools software package (this is not R). By default, the name is assumed to be <bamfile>.bai.}
\item{bedfile}{A file with aligned reads in .bed format.}
\item{bedGraphfile}{A file with signal values in .bedGraph format.}
\item{chrom.length.file}{
A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.}
\item{outputfolder}{Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.}
\item{binsizes}{Integer value(s) which will be used for the binning. If more than one value is given, output files will be produced for each bin size.}
\item{chromosomes}{If only a subset of the chromosomes should be binned, specify them here.}
\item{separate.chroms}{
If set to 'TRUE', separate output files will be written for each chromosome. The chromosome name will be included in the filename. If set to 'FALSE', only a single output file containing all the chromosomes will be produced.
}
\item{save.as.RData}{
If set to 'FALSE', no output file will be written. Instead, a data.frame containing the binned data will be returned. Only the first binsize and the first (separate) chromosome will be processed in this case!
}
}

\details{
These functions convert aligned reads or signal values from various file formats into read counts in equidistant windows (bins). bam2binned() and bed2binned() use 'GenomicRanges::countOverlaps()' to calculate the read counts. Therefore, with small binsizes and large read lengths, one read fragment will often overlap more than one bin.
bedGraph2binned() sets the maximum signal value in a bin as value for that bin.
}
\value{
Return value is a data.frame with 4 columns. Columns 1 to 3 are chromosomal coordinates and column 4 is the read count in the corresponding bin. This data.frame will be either written to file or given as return value, depending on the option 'save.as.RData'. It will look like the following:
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
\keyword{bam}
\keyword{bed}
\keyword{bedGraph}
\keyword{binned}
\keyword{bin}
\keyword{window}
\references{}
\author{Aaron Taudt}
\note{}
\seealso{\code{\link{binned2wiggle}}}
\examples{
}
