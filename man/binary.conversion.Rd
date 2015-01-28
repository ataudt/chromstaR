\name{Decimal to binary conversion}
\alias{dec2bin}
\alias{bin2dec}

\title{
Conversion of decimal and binary states
}
\description{
Convert combinatorial states in decimal representation to combinatorial states in binary representation and vice versa.
}
\usage{
dec2bin(dec, ndigits=NULL, names=NULL)
bin2dec(bin)
}
\arguments{
\item{dec}{An integer vector.}
\item{ndigits}{The number of digits that the binary representation should have. If unspecified, the shortest possible representation will be chosen.}
\item{names}{A vector of characters. If specified, a data.frame with names will be returned instead of the binary matrix.}
\item{bin}{A matrix with only 0 and 1 (or TRUE and FALSE) as entries. One combinatorial state per row.}
}

\details{
\pkg{\link{chromstaR}} uses decimal numbers to represent combinatorial states of peaks. These functions serve as a convenient way to get from the efficient decimal representation to a more human-readable binary representation.
}
\value{
A vector of integers for \code{bin2dec} and a matrix of logicals with one state per row for \code{dec2bin}.
}
\keyword{binary}
\keyword{decimal}
\keyword{conversion}
\author{Aaron Taudt}
\examples{
}

