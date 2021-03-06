% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binReads.R
\name{binReads}
\alias{binReads}
\alias{binning}
\title{Convert aligned reads from various file formats into read counts in equidistant bins}
\usage{
binReads(
  file,
  experiment.table = NULL,
  ID = NULL,
  assembly,
  bamindex = file,
  chromosomes = NULL,
  pairedEndReads = FALSE,
  min.mapq = 10,
  remove.duplicate.reads = TRUE,
  max.fragment.width = 1000,
  blacklist = NULL,
  binsizes = 1000,
  stepsizes = binsizes/2,
  reads.per.bin = NULL,
  bins = NULL,
  variable.width.reference = NULL,
  use.bamsignals = TRUE,
  format = NULL
)
}
\arguments{
\item{file}{A file with aligned reads. Alternatively a \code{\link{GRanges-class}} with aligned reads.}

\item{experiment.table}{An \code{\link{experiment.table}} containing the supplied \code{file}. This is necessary to uniquely identify the file in later steps of the workflow. Set to \code{NULL} if you don't have it (not recommended).}

\item{ID}{Optional ID to select a row from the \code{experiment.table}. Only necessary if the experiment table contains the same file in multiple positions in column 'file'.}

\item{assembly}{Please see \code{\link[GenomeInfoDb]{getChromInfoFromUCSC}} for available assemblies. Only necessary when importing BED files. BAM files are handled automatically. Alternatively a data.frame with columns 'chromosome' and 'length'.}

\item{bamindex}{BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.}

\item{chromosomes}{If only a subset of the chromosomes should be binned, specify them here.}

\item{pairedEndReads}{Set to \code{TRUE} if you have paired-end reads in your BAM files (not implemented for BED files).}

\item{min.mapq}{Minimum mapping quality when importing from BAM files. Set \code{min.mapq=0} to keep all reads.}

\item{remove.duplicate.reads}{A logical indicating whether or not duplicate reads should be removed.}

\item{max.fragment.width}{Maximum allowed fragment length. This is to filter out erroneously wrong fragments due to mapping errors of paired end reads.}

\item{blacklist}{A \code{\link{GRanges-class}} or a bed(.gz) file with blacklisted regions. Reads falling into those regions will be discarded.}

\item{binsizes}{An integer vector specifying the bin sizes to use.}

\item{stepsizes}{An integer vector specifying the step size. One number can be given for each element in \code{binsizes}, \code{reads.per.bin} and \code{bins} (in that order).}

\item{reads.per.bin}{Approximate number of desired reads per bin. The bin size will be selected accordingly.}

\item{bins}{A \code{\link{GRanges-class}} or a named \code{list()} with \code{\link{GRanges-class}} containing precalculated bins produced by \code{\link{fixedWidthBins}} or \code{\link{variableWidthBins}}. Names of the list must correspond to the binsize. If the list is unnamed, an attempt is made to automatically determine the binsize.}

\item{variable.width.reference}{A BAM file that is used as reference to produce variable width bins. See \code{\link{variableWidthBins}} for details.}

\item{use.bamsignals}{If \code{TRUE} the \pkg{\link[bamsignals]{bamsignals}} package is used for parsing of BAM files. This gives tremendous speed advantage for only one binsize but linearly increases for multiple binsizes, while \code{use.bamsignals=FALSE} has a binsize dependent runtime and might be faster if many binsizes are calculated.}

\item{format}{One of \code{c('bed','bam','GRanges',NULL)}. With \code{NULL} the format is determined automatically from the file ending.}
}
\value{
If only one bin size was specified for option \code{binsizes}, the function returns a single \code{\link{GRanges-class}} object with meta data column 'counts' that contains the read count. If multiple \code{binsizes} were specified , the function returns a named \code{list()} of \link{GRanges-class} objects.
}
\description{
Convert aligned reads in .bam or .bed(.gz) format into read counts in equidistant windows.
}
\details{
Convert aligned reads from .bam or .bed(.gz) files into read counts in equidistant windows (bins). This function uses \code{GenomicRanges::countOverlaps} to calculate the read counts, or alternatively \code{bamsignals::bamProfile} if option \code{use.bamsignals} is set (only effective for .bam files).
}
\examples{
## Get an example BAM file with ChIP-seq reads
file <- system.file("extdata", "euratrans",
                   "lv-H3K27me3-BN-male-bio2-tech1.bam",
                    package="chromstaRData")
## Bin the file into bin size 1000bp
data(rn4_chrominfo)
data(experiment_table)
binned <- binReads(file, experiment.table=experiment_table,
                  assembly=rn4_chrominfo, binsizes=1000,
                  stepsizes=500, chromosomes='chr12')
print(binned)

}
