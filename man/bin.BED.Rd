\name{bin.BED}
\alias{bin.BED}
\title{
Bin a BED file
}
\description{
The function will bin the reads of an input BED file into a specified binsize. It takes a text file in BED format as input and produces output files with the binned data. The output files will be named as the input file with additional information (binsize, chromosome) appended to the filename.
}
\usage{
bin.BED(bedfile, reference.genome.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE)
}
\arguments{
\item{bedfile}{
A file with mapped reads in BED format. Only the first 3 columns with the genomic coordinates will be used. The file could look like this for example:
\tabular{llllll}{
chr1 \tab 28560 \tab 28596 \tab N \tab 1000 \tab + \cr
chr1 \tab 28563 \tab 28599 \tab N \tab 1000 \tab - \cr
chr1 \tab 28571 \tab 28607 \tab N \tab 1000 \tab + \cr
chr1 \tab 61828 \tab 61864 \tab N \tab 1000 \tab - \cr
chr1 \tab 102946 \tab 102982 \tab N \tab 1000 \tab - \cr
chr1 \tab 235666 \tab 235702 \tab N \tab 1000 \tab - \cr
... \tab ... \tab ... \tab ... \tab ... \tab ... \cr
}
}
\item{reference.genome.file}{
A file which contains the chromosome lengths (in bp). The first column contains the chromosome name and the second column the length.
\tabular{ll}{
chr1 \tab 249250621 \cr
chr2 \tab 243199373 \cr
chr3 \tab 198022430 \cr
chr4 \tab 191154276 \cr
... \tab ... \cr
}
}
\item{outputfolder}{
Folder to which the binned data will be saved. If the specified folder does not exist, it will be created.
}
\item{binsizes}{
Integer value(s) which will be used for the binning.
}
\item{chromosomes}{
Chromosome name(s) which will be used for the binning.
}
\item{separate.chroms}{
If set to 'TRUE', separate output files will be written for each chromosome. The chromosome name will be included in the filename. If set to 'FALSE', only a single output file containing all the chromosomes will be produced.
}
\item{save.as.RData}{
If set to 'FALSE', no output file will be written. Instead, a data.frame containing the binned data will be returned. Only the first binsize and the first (separate) chromosome will be processed in this case!
}
}

\details{
The function takes a text file in BED format (\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}) as input and converts the mapped reads into binned representation. The following two tables should make clear what the function does:

The input file is in BED format and only the first three columns will be considered:
\tabular{llllll}{
chr1 \tab 28560 \tab 28596 \tab N \tab 1000 \tab + \cr
chr1 \tab 28563 \tab 28599 \tab N \tab 1000 \tab - \cr
chr1 \tab 28571 \tab 28607 \tab N \tab 1000 \tab + \cr
chr1 \tab 61828 \tab 61864 \tab N \tab 1000 \tab - \cr
chr1 \tab 102946 \tab 102982 \tab N \tab 1000 \tab - \cr
chr1 \tab 235666 \tab 235702 \tab N \tab 1000 \tab - \cr
... \tab ... \tab ... \tab ... \tab ... \tab ... \cr
}
The binned data will then look like the following:
\tabular{rrrr}{
chromosome \tab start \tab end \tab reads \cr
... \tab ... \tab ... \tab ... \cr
chr1 \tab 28000 \tab 28199 \tab 0 \cr
chr1 \tab 28200 \tab 28399 \tab 0 \cr
chr1 \tab 28400 \tab 28599 \tab 3 \cr
chr1 \tab 28600 \tab 28799 \tab 0 \cr
... \tab ... \tab ... \tab ... \cr
}
}
\value{
A data.frame containing the binned data. This data.frame will be either written to file or given as return value, depending on the option 'save.as.RData'.
}
\keyword{BED}
\keyword{bin}
\references{\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}}
\author{Aaron Taudt}
\note{}
\seealso{}
\examples{
}
