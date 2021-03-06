% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCombinations.R
\name{getCombinations}
\alias{getCombinations}
\title{Get combinations}
\usage{
getCombinations(gr)
}
\arguments{
\item{gr}{A \code{\link[GenomicRanges]{GRanges-class}} object from which the meta-data columns containing combinations will be extracted.}
}
\value{
A DataFrame.
}
\description{
Get a DataFrame with combinations from a \code{\link[GenomicRanges]{GRanges-class}} object.
}
\examples{
### Get an example multiHMM ###
file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
                    package="chromstaR")
model <- get(load(file))
### Get the combinations
bin.combs <- getCombinations(model$bins)
print(bin.combs)
seg.combs <- getCombinations(model$segments)
print(seg.combs)

}
