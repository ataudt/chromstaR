#ifndef LOGHMM_H
#define LOGHMM_H

#include "utility.h"
#include "densities.h"
#include "logging.h"
#include <omp.h>

class LogHMM
{

	public:
		int T; ///< length of observed sequence
		int N; ///< number of states
		int Nmod; ///< number of modifications / marks
		double** A; ///< matrix [N x N] of transition probabilities
		double** logA; ///< logarithm of matrix [N x N] of transition probabilities
		double* logproba; ///< initial probabilities (length N)
		double logP; ///< loglikelihood
		int* O; ///< vector [T] of observations
		int** multiO; ///< matrix [Nmod x T] of observations
		vector<Density*> densityFunctions; ///< density functions for each state
		
		LogHMM();
		LogHMM(int T, int N); ///< constructor for the univariate HMM
		LogHMM(int T, int N, int Nmod); ///< constructor for the multivariate HMM
		~LogHMM();
		void initialize_transition_probs(double* initial_A, bool use_initial_params);
		void initialize_proba(double* initial_proba, bool use_initial_params);
		void get_posteriors(double** post);
		void baumWelch(int* maxiter, int* maxtime, double* eps);
		void check_for_state_swap();
		void calc_weights(double* weights);
// 		void viterbi(int* path, int recompute);
		double get_proba(int i);
        
	private:
		void forward();
		void backward();
		void calc_sumgamma();
		void calc_sumxi();
		void calc_loglikelihood();
		void computeDensities();
		void print_uni_params();
		void print_multi_params();
		void print_uni_iteration(int iteration);
		void print_multi_iteration(int iteration);
		double** logalpha; ///< matrix[T x N] of forward probabilities
		double** logbeta; ///<  matrix[T x N] of backward probabilities
		double** logdensities; ///< matrix[N x T] of precomputed logdensity values
		double* sumgamma; ///< vector[N] of sum of posteriors (gamma values)
		double** sumxi; ///< matrix[N x N] of xi values
		double** gamma; ///< matrix[N x T] of posteriors
		double dlogP; ///< difference in loglikelihood from one iteration to the next
		time_t baumWelchStartTime_sec; ///< start time of the Baum-Welch in sec
		int baumWelchTime_real; ///< elapsed time from start of the 0th iteration
		int sumdiff_state1; ///< sum of the difference in the state 1 assignments from one iteration to the next
		double sumdiff_posterior; ///< sum of the difference in posterior (gamma) values from one iteration to the next
// 		int* num_nonzero_A_into_state; ///< vector[N] number of non-zero transition into each state for exploiting transition matrix sparsity
// 		int** index_nonzero_A_into_state; ///< matrix[N x N] of indices of non-zero transitions into each state for exploiting transition matrix sparsity
// 		double transition_cutoff; ///< cutoff value for considering transition values as zero
// 		double sparsity_cutoff; ///< cutoff value for using sparsity exploit in forward() and backward() step. If <= 0, sparsity exploit is never used, if >1, sparsity exploit is always used.
		whichvariate xvariate;
};

#endif
