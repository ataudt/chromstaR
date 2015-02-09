\name{Convert aligned reads to read counts}
\alias{bam2binned}
\alias{bed2binned}
\alias{bedGraph2binned}
\alias{binning}

\title{
Convert aligned reads from various file formats into read counts in equidistant bins
}
\description{
Convert aligned reads in .bam or .bed format into read counts in equidistant windows. Convert signal values in .bedGraph format to signal counts in equidistant windows.
}
\usage{
bam2binned(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE)

bed2binned(bedfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE)

bedGraph2binned(bedGraphfile, assembly, chrom.length.file=NULL, outputfolder="binned_data", binsizes=500, chromosomes=NULL, separate.chroms=FALSE, save.as.RData=TRUE)
}
\arguments{
\item{bamfile}{A file with aligned reads in .bam format.}
\item{bamindex}{A bam-index file (.bai). If you miss this file, generate it with the command 'samtools index bamfile.bam' of the SAMtools software package (this is not R). By default, the name is assumed to be <bamfile>.bai.}
\item{bedfile}{A file with aligned reads in .bed format.}
\item{bedGraphfile}{A file with signal values in .bedGraph format.}
\item{assembly}{An assembly to specify the chromosome lengths. One of \code{c('hg19','hg18')}.}
\item{chrom.length.file}{
A file which contains the chromosome lengths in basepairs. The first column contains the chromosome name and the second column the length (see also \code{\link{chrom.length.file}}.}
\item{outputfolder}{Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.}
\item{binsizes}{Integer value(s) which will be used for the binning. If more than one value is given, output files will be produced for each bin size.}
\item{chromosomes}{If only a subset of the chromosomes should be binned, specify them here.}
\item{separate.chroms}{
If set to 'TRUE', separate output files will be written for each chromosome. The chromosome name will be included in the filename. If set to 'FALSE' (default), only a single output file containing all the chromosomes will be produced.
}
\item{save.as.RData}{
If set to 'FALSE', no output file will be written. Instead, a \link{GenomicRanges} object containing the binned data will be returned. Only the first binsize will be processed in this case and option \code{separate.chroms} will be set to 'FALSE').
}
}

\details{
These functions convert aligned reads or signal values from various file formats into read counts in equidistant windows (bins). bam2binned() and bed2binned() use 'GenomicRanges::countOverlaps()' to calculate the read counts. Therefore, with small binsizes and large read lengths, one read fragment will often overlap more than one bin.
bedGraph2binned() sets the maximum signal value in a bin as value for that bin.
}
\value{
The function produces a \link{GenomicRanges} object with one meta data column 'reads' that contains the read count. This binned data will be either written to file (\code{save.as.RData=TRUE}) or given as return value (\code{save.as.RData=FALSE}).
}
\keyword{bam}
\keyword{bed}
\keyword{bedGraph}
\keyword{binned}
\keyword{bin}
\keyword{window}
\author{Aaron Taudt}
\seealso{\link{chromstaR_binnedData}}
\examples{
## Get an example BED-file with ChIP-seq reads for H3K36me3 in brain tissue
bedfile <- system.file("extdata/brain/BI.Brain_Angular_Gyrus.H3K36me3.112.bed_chr22.gz", package="chromstaR")
## Bin the BED files into bin size 800bp
binned.data <- bed2binned(bedfile, assembly='hg19', binsize=800, save.as.RData=FALSE)
}
