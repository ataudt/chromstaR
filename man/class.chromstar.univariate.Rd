\name{chromstar univariate object}
\alias{chromstar.univariate.hmm}

\title{Univariate HMM object}
\description{
The univariate HMM object is output of the function \code{\link{univariate.from.binned.data}} and is basically a list with various entries. The class() attribute of this list was set to "chromstar.univariate.hmm". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
}
\value{
\item{reads}{Vector of read counts for each bin.}
\item{num.bins}{Number of bins. Same as \emph{length(hmm$reads)}.}
\item{num.states}{Number of states (3) that define this HMM. The states are 'zero-inflation', 'unmodified' and 'modified'.}
\item{num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \emph{eps}.}
\item{time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \emph{eps}.}
\item{loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \emph{eps})}
\item{posteriors}{Matrix of posterior probabilities for each bin (rows) and state (columns).}
\item{A}{Matrix of transition probabilities for each state into each state.}
\item{proba}{Probabilities for the first bin. Same as \emph{hmm$posteriors[1,]}.}
\item{loglik}{Final loglikelihood after the last iteration.}
\item{weights}{Weight for each component. Same as \emph{apply(hmm$posteriors,2,mean)}.}
\item{A.initial}{Initial \emph{A} at the beginning of the Baum-Welch.}
\item{proba.initial}{Initial \emph{proba} at the beginning of the Baum-Welch.}
\item{num.threads}{Number of threads that were used in the Baum-Welch.}
\item{error}{This is 0 unless an error occurred during the Baum-Welch.}
\item{eps}{Convergence threshold for the Baum-Welch.}
\item{distributions}{Estimated parameters of the Negative Binomials.}
\item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
\item{coordinates}{Genomic coordinates of the bins.}
\item{states}{Classification of the bins into 'unmodified' and 'modified'. If you also want the 'zero-inflation' as a separate state, you can call 'get.states(hmm, separate.zeroinflation=T) -> hmm$states'.}
}

\keyword{univariate}
\keyword{hmm}
\seealso{\code{\link{chromstar.multivariate.hmm}}}

