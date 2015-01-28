\name{stateBrewer}
\alias{stateBrewer}

\title{
Obtain combinatorial states from specification
}
\description{
This function returns all combinatorial (decimal) states that are consistent with a given abstract specification.
}
\usage{
stateBrewer(statespec, diffstatespec=NULL)
}
\arguments{
\item{statespec}{A vector composed of any combination of the following entries: \code{0, 1, 'x', 'r[...]'}, where [\dots] can be any string.
	\itemize{
		\item \code{0}: sample is 'unmodified'
		\item \code{1}: sample is 'modified'
		\item \code{'x'}: sample can be both 'unmodified' or 'modified'
		\item \code{'rA'}: all samples in group A have to be in the same state
		\item \code{'rB'}: all samples in group B have to be in the same state
		\item \code{'r[...]'}: all samples in group [\dots] have to be in the same state
	}
}
\item{diffstatespec}{A vector composed of any combination of the following entries: \code{'x', 'd[]'}, where [] can be any string.
	\itemize{
		\item \code{'x'}: sample can be both 'unmodified' or 'modified'
		\item \code{'dA'}: at least one sample in group A has to be different from the other samples in group A 
		\item \code{'dB'}: at least one sample in group B has to be different from the other samples in group B 
		\item \code{'d[...]'}: at least one sample in group [\dots] has to be different from the other samples in group [\dots] 
	}
}
}
\details{
The binary modification state (unmodified=0 or modified=1) of multiple ChIP-seq samples defines a (decimal) combinatorial state such as:
\tabular{ccccccc}{
 \tab sample1 \tab sample2 \tab sample3 \tab sample4 \tab sample5 \tab combinatorial state \cr
bin1 \tab 0 \tab 0 \tab 1 \tab 0 \tab 0 \tab 4 \cr
bin2 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
bin3 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \tab 10 \cr
bin4 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 15 \cr
bin5 \tab 0 \tab 0 \tab 1 \tab 0 \tab 1 \tab 5 \cr
}
}

\value{
A integer vector with (decimal) combinatorial states following the given specification.
}
\keyword{combinatorial}
\author{Aaron Taudt}
\examples{
# Get all combinatorial states where sample1=0, sample2=1, sample3=(0 or 1), sample4=sample5
stateBrewer(statespec=c(0,1,'x','rA','rA'))

# Get all combinatorial states where sample1=sample2=sample3, sample4=sample5
stateBrewer(statespec=c('rA','rA','rA','rB','rB'))

# Get all combinatorial states where sample1=sample5, sample2=sample3=1, sample4=(0 or 1)
stateBrewer(statespec=c('rA',1,1,'x','rA'))

# Get all combinatorial states where sample1 != sample2, sample3 != sample4, sample5=(0 or 1)
stateBrewer(statespec=c('x','x','x','x','x'), diffstatespec=c('dA','dA','dB','dB','x'))
# To check whether you chose the correct specification you can use
dec2bin(stateBrewer(statespec=c('x','x','x','x','x'), diffstatespec=c('dA','dA','dB','dB','x')), ndigits=5)
}
