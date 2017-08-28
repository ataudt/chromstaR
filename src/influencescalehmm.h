#ifndef INFLUENCESCALEHMM_H
#define INFLUENCESCALEHMM_H

#include <R.h> // R_CheckUserInterrupt()
#include <vector> // storing density functions
#include <time.h> // time(), difftime()
#include "utility.h"
#include "densities.h"
#include <string> // strcmp
#include <cstdlib> // std::abs

#ifdef _OPENMP
#include <omp.h>
#endif
// #if defined TARGET_OS_MAC || defined __APPLE__
// #include <libiomp/omp.h> // parallelization options on mac
// #elif defined __linux__ || defined _WIN32 || defined _WIN64
// #include <omp.h> // parallelization options
// #endif

class InfluenceScaleHMM  {

	public:
		// Constructor and Destructor
		InfluenceScaleHMM();
		//InfluenceScaleHMM(int T, int N,  int verbosity);
		InfluenceScaleHMM(int T, int N, int Nmod, int verbosity);
		~InfluenceScaleHMM();

		// Member variables
		std::vector<Density*> densityFunctions; ///< density functions for each state

		// Methods
		void initialize_transition_probs(double* initial_A, bool use_initial_params);
		void initialize_proba(double* initial_proba, bool use_initial_params);
		//void initialize_influence();
		void initialize_tiestrength(double* initial_tiestrenght, bool use_initial_params);
		void baumWelch(int* maxiter, int* maxtime, double* eps);
		//void check_for_state_swap();
		std::vector< std::vector< double > > calc_weights();
		void calc_weights(double* weights);

		// Getters and Setters
		int get_N();
		int get_T();
		void get_posteriors(double**** post);
		double get_posterior(int c1, int iN, int t);
		double get_density(int c1, int iN, int t);
		double get_proba(int c1, int i);
		double get_A(int c1, int c2, int i, int j);
		double get_logP();
		void set_cutoff(int cutoff);

	private:
		// Member variables
		int verbosity; ///< verbosity level
		int T; ///< length of observed sequence
		int N; ///< number of states
		int C;
		int cutoff; ///< a cutoff for observations
		int Nmod; ///< number of modifications / marks
		double**** A; ///< matrix [N x N] of transition probabilities
		double** proba; ///< initial probabilities (length N)
		double** weights;
		double logP; ///< loglikelihood
		int* O; ///< vector [T] of observations
		int** multiO; ///< matrix [Nmod x T] of observations
		double** scalefactoralpha; ///< vector[T] of scaling factors
		double*** scalealpha; ///< matrix [T x N] of forward probabilities
		double*** scalebeta; ///<  matrix [T x N] of backward probabilities
		double*** densities; ///< matrix [N x T] of density values
//	double** tdensities; ///< matrix [T x N] of density values, for use in multivariate !increases speed, but on cost of RAM usage and that seems to be limiting
		double** sumgamma; ///< vector[N] of sum of posteriors (gamma values)
		double**** influence;
		double** tiestrength;
		double**** sumxi; ///< matrix[N x N] of xi values
		double** tiesumxi;
		double* tiesumxitotal;
		double*** tie_onestatus_sumxi;
		double*** gamma; ///< matrix[N x T] of posteriors
		bool* states_prev; ///< vector[T] with modification state of last state ('modified') of previous iteration
		double dlogP; ///< difference in loglikelihood from one iteration to the next
		time_t baumWelchStartTime_sec; ///< start time of the Baum-Welch in sec
		int baumWelchTime_real; ///< elapsed time from start of the 0th iteration
		int sumdiff_state_last; ///< sum of the difference in the state 1 assignments from one iteration to the next (only univariate case)
		whichvariate xvariate; ///< enum which stores if UNIVARIATE or MULTIVARIATE

		// Methods
		void calc_influence();
		void calc_tiestrength();
		void forward();
		void backward();
		void calc_sumgamma();
		void calc_sumxi();
		void calc_tiesumxi();
		void calc_loglikelihood();
		void calc_densities();
		void print_uni_iteration(int iteration);
		void print_uni_params();
		void print_multi_iteration(int iteration);
		void print_multi_params();
};

#endif
