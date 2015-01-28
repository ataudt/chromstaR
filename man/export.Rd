\name{Export to Genome Browser}
\alias{export.unihmm2bed}
\alias{export.unihmm2wiggle}
\alias{export.multihmm2bed}
\alias{export.multihmm2wiggle}
\alias{export.GRanges2bed}
\alias{export.binned2bedGraph}
\alias{export.binned2bed}
\alias{export.binned2wiggle}
\alias{export}

\title{
Export various \pkg{\link{chromstaR}} objects into genome browser viewable files.
}
\description{
Export peak calls or read counts as genome browser viewable files.
}
\usage{
export.unihmm2bed(hmm.list, only.modified=TRUE, filename="view_me_in_genome_browser")

export.unihmm2wiggle(hmm.list, filename="view_me_in_genome_browser")

export.multihmm2bed(multi.hmm, separate.tracks=TRUE, exclude.states=0, include.states=NULL, filename="view_me_in_genome_browser")

export.multihmm2wiggle(multi.hmm, filename="view_me_in_genome_browser")

export.GRanges2bed(gr, filename="view_me_in_genome_browser", separate.tracks=FALSE)

export.binned2bedGraph(binned.data.list, filename="view_me_in_genome_browser")

export.binned2wiggle(binned.data.list, filename="view_me_in_genome_browser")
}
\arguments{
\item{hmm.list}{A list of \code{\link{chromstaR_univariateHMM}} objects, e.g. \code{list(hmm1, hmm2, ...)}.}
\item{multi.hmm}{A \code{\link{chromstaR_multivariateHMM}} object.}
\item{gr}{A \link{GRanges} object.}
\item{binned.data.list}{A list of \link{GRanges} objects containing binned read counts, e.g. \code{list(binned1, binned2, ...)}.}
\item{filename}{Name of the output file. The appropriate ending will be appended.}
\item{only.modified}{If set to TRUE, only coordinates of regions in state 'modified' (=peaks) will be exported.}
\item{separate.tracks}{Export the different combinatorial states in one track with different colors (FALSE) or in separate tracks (TRUE).}
\item{exclude.states}{An integer vector with combinatorial states that shall not be exported.}
\item{include.states}{An integer vector with combinatorial states that shall be exported. If not specified, all available states except the excluded ones will be exported.}
}

\details{
These functions export various \pkg{\link{chromstaR}} objects like \link{chromstaR_binnedData}, \code{\link{chromstaR_univariateHMM}}s and \code{\link{chromstaR_multivariateHMM}}s into files which can be uploaded into a Genome Browser. Choose the wiggle-format if you want to export binned read counts from your object. Choose bed-format if you want to export peak calls.
}
\value{
A file which can be uploaded into a genome browser.
}
\keyword{export}
\keyword{genome}
\keyword{browser}
\keyword{wiggle}
\keyword{bed}
\author{Aaron Taudt}
\examples{
}
