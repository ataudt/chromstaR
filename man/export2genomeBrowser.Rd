\name{Export to Genome Browser}
\alias{univariate2bed}
\alias{binned2wiggle}
\alias{multivariate2bed}

\title{
Export various \code{\link{chromstar}} objects into Genome Browser viewable files.
}
\description{
Export binned read counts or the results of the Hidden Markov Models as Genome Browser viewable files.
}
\usage{
binned2wiggle(binned.data.list, file="view_me_in_genome_browser")

univariate2bed(uni.hmm.list, file="view_me_in_genome_browser", threshold=0.5, separate.zeroinflation=FALSE, chrom.length.file=NULL)

multivariate2bed(multi.hmm, separate.tracks=FALSE, exclude.state.zero=TRUE, numstates=NULL, file="view_me_in_genome_browser", chrom.length.file=NULL)
}
\arguments{
\item{binned.data.list}{A named list of data.frames with binned read counts.}
\item{file}{Name of the output file.}
\item{uni.hmm.list}{A list of univariate HMMs, produced by the function \code{\link{univariate.from.binned.data}}.}
\item{threshold}{
Threshold for calling a bin 'modified' or 'unmodified'. If the posterior for state 'modified' is greater than the threshold, then this bin will be called 'modified'. The default of 0.5 should be best in most cases.
}
\item{separate.zeroinflation}{
If set to FALSE, the 'zero-inflation' and 'unmodified' state will be merged into 'unmodified'.
}
\item{chrom.length.file}{
An optional file with chromosome lengths. If this file is given, the chromosome boundaries of the input data are checked and corrected in case they are wrong (see \code{\link{chrom.length.file}} for a description of this file).
}
\item{multi.hmm}{A multivariate HMM, produced by the function \code{\link{multivariate.from.univariate.results}}.}
\item{separate.tracks}{If set to TRUE, the combinatorial states will be displayed as separate tracks instead of a single one.}
\item{exclude.state.zero}{If set to TRUE, combinatorial state 0 (everything unmodified) is excluded from the export.}
\item{numstates}{Number of combinatorial states to export (order as in 'multi.hmm$comb.states')}
}

\details{
These functions export various \code{\link{chromstar}} objects like \code{\link{binned data}}, \code{\link{univariate HMM results}} and \code{\link{multivariate HMM results}} into files which can be uploaded into a Genome Browser. The output file formats (wiggle or bed) depend on the object to be exported and are chosen to represent the data most efficiently.
}
\value{
A file which can be uploaded into a Genome Browser.
}
\keyword{export}
\keyword{genome}
\keyword{browser}
\keyword{wiggle}
\keyword{bed}
\references{}
\author{Aaron Taudt}
\note{}
\seealso{}
\examples{
}
