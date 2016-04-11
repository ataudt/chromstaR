#' Binned read counts
#'
#' A \link{GRanges} object which contains binned read counts as meta data column \code{counts}. It is output of the \code{\link{binReads}} function.
#' @name binned.data
NULL

#' Univariate HMM object
#'
#' The univariate HMM object is output of the function \code{\link{callPeaksUnivariate}} and is basically a list with various entries. The class() attribute of this list was set to "uniHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
#'
#' @return
#' A \code{list()} with the following entries:
#' \item{ID}{An identifier that is used in various \pkg{\link{chromstaR}} functions.}
#' \item{bins}{
#' A \link{GRanges} object containing the genomic bin coordinates, their read count, (optional) posteriors and state classification.
#' }
#' \item{segments}{
#' A \link{GRanges} object containing regions and their state classification.
#' }
#' \item{weights}{Weight for each component. Same as \code{apply(hmm$posteriors,2,mean)}.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin. Same as \code{hmm$posteriors[1,]}.}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Estimated parameters of the emission distributions.}
#' \item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
#' \item{FDR}{False discovery rate (default=0.5).}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#' @seealso \code{\link{callPeaksUnivariate}}, \code{\link{multiHMM}}
#' @name uniHMM
#' @aliases uni.hmm
NULL

#' Multivariate HMM object
#'
#' The multivariate HMM object is output of the function \code{\link{callPeaksMultivariate}} and is basically a list with various entries. The class() attribute of this list was set to "multiHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' and 'hmm$'.
#' 
#' @return
#' A \code{list()} with the following entries:
#' \item{IDs}{IDs of the input univariate HMMs.}
#' \item{bins}{
#' A \link{GRanges} object containing the genomic bin coordinates, their read count, (optional) posteriors and state classification.
#' }
#' \item{segments}{
#' A \link{GRanges} object containing regions and their state classification.
#' }
#' \item{weights}{Weight for each component. Same as \code{apply(hmm$posteriors,2,mean)}.}
#' \item{weights.univariate}{Weights of the univariate HMMs.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin. Same as \code{hmm$posteriors[1,]}.}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Emission distributions used for this model.}
#' \item{FDR}{False discovery rate. NULL means that the state with maximum posterior probability was chosen, irrespective of its absolute probability (default=NULL).}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{correlation.matrix}{Correlation matrix of transformed reads.}
#' @seealso \code{\link{callPeaksMultivariate}}, \code{\link{uniHMM}}
#' @name multiHMM
#' @aliases multi.hmm
NULL
