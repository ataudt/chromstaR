#' Experiment data table
#'
#' A \code{data.frame} specifying the structure of the experiment.
#'
#' @format A \code{data.frame} with columns 'file', 'mark', 'condition', 'replicate' and 'pairedEndReads'. Avoid the use of special characters like '-' or '+' as this will confuse the internal file management.
#' @name experiment.table
#' @examples
#'data(experiment_table)
#'print(experiment_table)
#'
NULL


#' Multivariate HMM for demonstration purposes
#'
#' A \code{\link{multiHMM}} object for demonstration purposes in examples of package \pkg{\link{chromstaR}}.
#'
#' @docType data
#' @name multivariate_model
#' @format A \code{\link{multiHMM}} object.
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
NULL


#' Combined multivariate HMM for demonstration purposes
#'
#' A \code{\link{combinedMultiHMM}} object for demonstration purposes in examples of package \pkg{\link{chromstaR}}.
#'
#' @docType data
#' @name combined_model
#' @format A \code{\link{combinedMultiHMM}} object.
#' @examples
#'## Get an example combinedMultiHMM
#'file <- system.file("data","combined_mode-differential.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
NULL


#' Gene coordinates for rn4
#'
#' A data.frame containing gene coordinates and biotypes of the rn4 assembly.
#'
#' @docType data
#' @name genes_rn4
#' @format A data.frame.
#' @examples
#'data(genes_rn4)
#'head(genes_rn4)
NULL


#' chromstaR objects
#'
#' @description
#' \pkg{\link{chromstaR}} defines several objects.
#' \itemize{
#' \item \code{\link{uniHMM}}: Returned by \code{\link{callPeaksUnivariate}}.
#' \item \code{\link{multiHMM}}: Returned by \code{\link{callPeaksMultivariate}} and \code{\link{callPeaksReplicates}}.
#' \item \code{\link{combinedMultiHMM}}: Returned by \code{\link{combineMultivariates}}.
#' }
#'
#' @name chromstaR-objects
NULL


#' Binned read counts
#'
#' A \code{\link[GenomicRanges]{GRanges}} object which contains binned read counts as meta data column \code{counts}. It is output of the \code{\link{binReads}} function.
#' @name binned.data
NULL


#' Univariate HMM object
#'
#' The univariate HMM object is output of the function \code{\link{callPeaksUnivariate}} and is a \code{list()} with various entries. The class() attribute of this list was set to "uniHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' or 'hmm$'.
#'
#' @return
#' A \code{list()} with the following entries:
#' \item{ID}{An identifier that is used in various \pkg{\link{chromstaR}} functions.}
#' \item{bins}{A \code{\link[GenomicRanges]{GRanges}} object containing the genomic bin coordinates, their read count, (optional) posteriors and state classification.}
#' \item{segments}{Same as \code{bins}, but consecutive bins with the same state are collapsed into segments.}
#' \item{weights}{Weight for each component. Same as \code{apply(hmm$posteriors,2,mean)}.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin. Same as \code{hmm$posteriors[1,]}.}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Estimated parameters of the emission distributions.}
#' \item{distributions.initial}{Distribution parameters at the beginning of the Baum-Welch.}
#' \item{post.cutoff}{Cutoff for posterior probabilities to call peaks (default=0.5).}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$max.mean}{Value of parameter \code{max.mean}.}
#' \item{convergenceInfo$read.cutoff}{Cutoff value for read counts.}
#' @seealso \code{\link{callPeaksUnivariate}}, \code{\link{multiHMM}}, \code{\link{combinedMultiHMM}}
#' @name uniHMM
#' @aliases uni.hmm
NULL

#' Multivariate HMM object
#'
#' The multivariate HMM object is output of the function \code{\link{callPeaksMultivariate}} and is a \code{list()} with various entries. The class() attribute of this list was set to "multiHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' or 'hmm$'.
#' 
#' @return
#' A \code{list()} with the following entries:
#' \item{IDs}{IDs of the input univariate HMMs.}
#' \item{bins}{A \code{\link[GenomicRanges]{GRanges}} object containing the genomic bin coordinates, their read count, (optional) posteriors and state classification.}
#' \item{segments}{Same as \code{bins}, but consecutive bins with the same state are collapsed into segments.}
#' \item{mapping}{A named vector giving the mapping from decimal combinatorial states to human readable combinations.}
#' \item{weights}{Weight for each component. Same as \code{apply(hmm$posteriors,2,mean)}.}
#' \item{weights.univariate}{Weights of the univariate HMMs.}
#' \item{transitionProbs}{Matrix of transition probabilities from each state (row) into each state (column).}
#' \item{transitionProbs.initial}{Initial \code{transitionProbs} at the beginning of the Baum-Welch.}
#' \item{startProbs}{Probabilities for the first bin. Same as \code{hmm$posteriors[1,]}.}
#' \item{startProbs.initial}{Initial \code{startProbs} at the beginning of the Baum-Welch.}
#' \item{distributions}{Emission distributions used for this model.}
#' \item{post.cutoff}{False discovery rate. NULL means that the state with maximum posterior probability was chosen, irrespective of its absolute probability (default=NULL).}
#' \item{convergenceInfo}{Contains information about the convergence of the Baum-Welch algorithm.}
#' \item{convergenceInfo$eps}{Convergence threshold for the Baum-Welch.}
#' \item{convergenceInfo$loglik}{Final loglikelihood after the last iteration.}
#' \item{convergenceInfo$loglik.delta}{Change in loglikelihood after the last iteration (should be smaller than \code{eps})}
#' \item{convergenceInfo$num.iterations}{Number of iterations that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{convergenceInfo$time.sec}{Time in seconds that the Baum-Welch needed to converge to the desired \code{eps}.}
#' \item{correlation.matrix}{Correlation matrix of transformed reads.}
#' @seealso \code{\link{callPeaksMultivariate}}, \code{\link{uniHMM}}, \code{\link{combinedMultiHMM}}
#' @name multiHMM
#' @aliases multi.hmm
#' @examples
#'## Get an example multiHMM
#'file <- system.file("data","multivariate_mode-combinatorial_condition-SHR.RData",
#'                     package="chromstaR")
#'model <- get(load(file))
NULL


#' Combined multivariate HMM object
#'
#' The multivariate HMM object is output of the function \code{\link{combineMultivariates}} and is a \code{list()} with various entries. The class() attribute of this list was set to "combinedMultiHMM". For a given hmm, the entries can be accessed with the list operators 'hmm[[]]' or 'hmm$'.
#' 
#' @return
#' A \code{list()} with the following entries:
#' \item{bins}{A \code{\link[GenomicRanges]{GRanges}} object containing genomic bin coordinates and human readable combinations for the combined \code{\link{multiHMM}} objects.}
#' \item{segments}{Same as \code{bins}, but consecutive bins with the same state are collapsed into segments.}
#' \item{segments.per.condition}{A \code{list} with segments for each condition separately.}
#' @seealso \code{\link{combineMultivariates}}, \code{\link{uniHMM}}, \code{\link{multiHMM}}
#' @name combinedMultiHMM
#' @aliases combinedHMM
NULL
