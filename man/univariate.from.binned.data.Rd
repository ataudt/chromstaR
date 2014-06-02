\name{univariate.from.binned.data}
\alias{univariate.from.binned.data}
\title{
Run an HMM for a single histone mark
}
\description{
Run an HMM for a single histone mark to determine regions which are either 'unmodified' or 'modified'. Input is a data.frame with \code{\link{binned.data}} and output is an \code{\link{chromstar.univariate.hmm}} list object.
}
\usage{
univariate.from.binned.data(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE, control=FALSE)
}
\arguments{
\item{binned.data}{A data.frame containing genomic coordinates and read counts (see \code{\link{binned.data}} for a more detailed description).}
\item{eps}{Convergence threshold for the Baum-Welch algorithm.}
\item{max.time}{The maximum running time for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.}
\item{max.it}{The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.}
\item{num.trials}{The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output. TODO: starting values are not random at the moment}
\item{eps.try}{If option num.trials is set to greater than 1, \option{eps.try} is used for the trial runs. If unset, \option{eps} is used.}
\item{num.threads}{Number of threads to use. Setting this to >1 may give increased performance.}
\item{output.if.not.converged}{
If set to FALSE, no output will be given in the case that the Baum-Welch did not converge to the desired \option{eps}. If set to TRUE, a model will be returned even if the Baum-Welch did not converge to the desired \option{eps}.
}
\item{filter.reads}{
If set to TRUE, then read counts which are above the 99.99\% quantile of the data are set to this value. A warning is issued in this case. This filtering can increase the speed of the Baum-Welch and is unlikely to affect the results.
}
\item{control}{
If set to TRUE, the binned data will be treated as control experiment. That means only state 'zero-inflation' and 'unmodified' will be used in the HMM.
}
}

\details{
The Hidden Markov Model which is used to classify the bins uses 3 states: state 'zero-inflation' with a delta function as emission densitiy (only zero read counts), 'unmodified' and 'modified' with Negative Binomials as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method.
TODO: insert paper
The function expects a data.frame with 4 columns (see arguments). Output of the function is a list object which contains the posteriors of the Baum-Welch likelihood maximization and the parameters for the Negative Binomial distributions. It also contains miscellaneous other attributes. All attributes can be accessed with the '$' or '[[ ]]' operator. To see which attributes are available, type "names(output)" (without quotation marks).
}
\value{
The output is a \code{\link{chromstar.univariate.hmm}} list object.
}
\keyword{Hidden Markov Model}
\keyword{Baum-Welch}
\keyword{univariate}
\references{}
\author{Aaron Taudt and Maria Colome Tatche}
\note{}
\seealso{}
\examples{
# Load the example dataset with genome coordinates and number of aligned reads for each bin
data(binned.data.example.1)
# Do the univariate HMM
model = univariate.from.binned.data(binned.data.example.1, eps=1, max.time=60)
# Look at the output
names(model)
# Select the posteriors
head(model$posteriors)
}
