\name{chromstar multivariate object}
\alias{chromstar.multivariate.hmm}
\alias{multi.hmm}

\title{Multivariate HMM object}
\description{
The multivariate HMM object is output of the function \code{\link{multivariate.from.univariate.results}} and is basically a list with various entries. The class() attribute of this list was set to "chromstar.multivariate.hmm". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
}
\value{
\item{reads}{Matrix of read counts for each bin (rows) and modification (columns).}
\item{num.bins}{Number of bins. Same as \option{nrow(hmm$reads)}.}
\item{num.states}{Number of combinatorial states that define this HMM.}
\item{num.modifications}{Number of modifications. Same as \option{ncol(hmm$reads)}.}
\item{comb.states}{Combinatorial states of the HMM.}
\item{correlation.matrix.inverse}{Inverse of the correlation matrix used for computing the densities.}
\item{determinant}{Determinant used for computing the densities.}
\item{num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \option{eps}.}
\item{time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \option{eps}.}
\item{loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \option{epsilon})}
\item{posteriors}{Matrix of posterior probabilities for each bin (rows) and state (columns).}
\item{A}{Matrix of transition probabilities for each state into each state.}
\item{proba}{Probabilities for the first bin. Same as \option{hmm$posteriors[1,]}.}
\item{loglik}{Final loglikelihood after the last iteration.}
\item{A.initial}{Initial \option{A} at the beginning of the Baum-Welch.}
\item{proba.initial}{Initial \option{proba} at the beginning of the Baum-Welch.}
\item{num.threads}{Number of threads that were used in the Baum-Welch.}
\item{error}{This is 0 unless an error occurred during the Baum-Welch.}
\item{coordinates}{Genomic coordinates of the bins.}
\item{eps}{Convergence threshold for the Baum-Welch.}
\item{distributions.univariate}{Parameters of the univariate Negative Binomials.}
\item{weights.univariate}{Weights of the univariate HMMs.}
\item{states}{Classification of the bins into combinatorial states.}
\item{correlation.matrix}{Correlation matrix of transformed reads.}
}

\keyword{multivariate}
\keyword{hmm}
\seealso{\code{\link{hmm2GRanges}}, \code{\link{multivariate2bed}}, \code{\link{chromstar.univariate.hmm}}}


