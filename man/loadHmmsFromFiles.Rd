\name{Load HMMs from files}
\alias{loadHmmsFromFiles}
\alias{loadMultiHmmsFromFiles}
\title{
Conveniently load Hidden Markov Models
}
\description{
These convenience functions save you the work of loading multiple saved \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} objects by hand. Input is a list of filenames with HMMs that you want to load, output is a list with the loaded HMM objects.
}
\usage{
loadHmmsFromFiles(filenames)
loadMultiHmmsFromFiles(filenames)
}
\arguments{
\item{filenames}{A character vector with files that contain \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} objects, respectively.}
}
\value{
A list with \code{\link{chromstaR_univariateHMM}} or \code{\link{chromstaR_multivariateHMM}} objects, respectively.
}
\keyword{load}
\author{Aaron Taudt}
\examples{
## Not run:
files.with.hmms <- list.files(folder.with.hmms, full=TRUE, pattern='binsize_800')
uni.HMMs <- loadHmmsFromFiles(files.with.hmms)
## End(Not run)
}
