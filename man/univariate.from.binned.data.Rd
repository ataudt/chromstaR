\name{univariate.from.binned.data}
\alias{univariate.from.binned.data}
\title{
Run an HMM for a single histone mark
}
\description{
Run an HMM for a single histone mark to determine regions which are either 'unmodified' or 'modified'. Input is a data.frame with \code{\link{binned.data}} and output is an \code{\link{uni.hmm}} list object.
}
\usage{
univariate.from.binned.data(binned.data, eps=0.001, max.time=-1, max.it=-1, num.trials=1, eps.try=NULL, num.threads=1, output.if.not.converged=FALSE, filter.reads=TRUE)
}
\arguments{
\item{binned.data}{A data.frame containing genomic coordinates and read counts (see \code{\link{binned.data}} for a more detailed description).}
\item{eps}{Convergence threshold for the Baum-Welch algorithm.}
\item{max.time}{The maximum running time for the Baum-Welch algorithm. If this time is reached, the Baum-Welch will terminate after the current iteration finishes. The default -1 is no limit.}
\item{max.it}{The maximum number of iterations for the Baum-Welch algorithm. The default -1 is no limit.}
\item{num.trials}{The number of trials to run the HMM. Each time, the HMM is seeded with different random initial values. The HMM with the best likelihood is given as output. TODO: starting values are not random at the moment}
\item{eps.try}{If option num.trials is set to greater than 1, \emph{eps.try} is used for the trial runs. If unset, \emph{eps} is used.}
\item{num.threads}{Number of threads to use. Setting this to >1 may give increased performance.}
\item{output.if.not.converged}{
If set to FALSE, no output will be given in the case that the Baum-Welch did not converge to the desired \emph{eps}. If set to TRUE, a model will be returned even if the Baum-Welch did not converge to the desired \emph{eps}.
}
\item{filter.reads}{
If set to TRUE, then read counts which are higher than 10 times the variance of the data are set to this value. A warning is issued in this case. This filtering is necessary because if the read counts are too high, the density becomes (numerically) zero for all states and the HMM crashes. It should not affect the results, though.
}
}

\details{
The Hidden Markov Model which is used to classify the bins uses 3 states: state 'zero-inflation' with a delta function as emission densitiy (only zero read counts), 'unmodified' and 'modified' with Negative Binomials as emission densities. A Baum-Welch algorithm is employed to estimate the parameters of the distributions. See our paper for a detailed description of the method.
TODO: insert paper
The function expects a data.frame with 4 columns (see arguments). Output of the function is a list object which contains the posteriors of the Baum-Welch likelihood maximization and the parameters for the Negative Binomial distributions. It also contains miscellaneous other attributes. All attributes can be accessed with the '$' or '[[ ]]' operator. To see which attributes are available, type "names(output)" (without quotation marks).
}
\value{
The output is an \code{\link{uni.hmm}} list object.
TODO: move these entries to the uni.hmm man page
\item{coordinates}{A data.frame with genomic coordinates. These are the same as in the input.}
\item{reads}{Number of aligned reads in the bins. Same as in input.}
\item{posteriors}{Posterior probabilities of the Baum-Welch estimation.}
\item{loglik}{Logarithm of the likelihood estimator (goodness of the fit).}
\item{iteration}{Number of iterations performed by the Baum-Welch.}
\item{time.in.sec}{Time in seconds spent in the Baum-Welch.}
\item{delta.loglik}{Change of the log(likelihood) in the last step. This should be smaller than the 'epsilon' in case of convergence.}
\item{epsilon}{Convergence threshold used for the Baum-Welch.}
\item{proba}{Estimated posteriors of the first bin.}
\item{A}{Estimated transition matrix.}
\item{distributions}{Estimated distribution parameters.}
\item{softweights}{Weights of the states. 'Soft' means that only posteriors where used to calculate the weights, as opposed to 'hard' assignments, where each bin is assigned to a state by using a threshold.}
\item{proba.initial}{Initial probabilities of the first bin.}
\item{A.initial}{Initial transition matrix.}
\item{distributions.initial}{Initial distribution parameters.}
\item{num.threads}{Number of threads used for the computations.}
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
