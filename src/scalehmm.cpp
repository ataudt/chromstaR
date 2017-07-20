#include "scalehmm.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------
ScaleHMM::ScaleHMM(int T, int N, int verbosity)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	//FILE_LOG(logDEBUG2) << "Initializing univariate ScaleHMM";
	this->xvariate = UNIVARIATE;
	this->verbosity = verbosity;
	this->T = T;
	this->N = N;
	this->A = CallocDoubleMatrix(N, N);
	this->scalefactoralpha = (double*) Calloc(T, double);
	this->scalealpha = CallocDoubleMatrix(T, N);
	this->scalebeta = CallocDoubleMatrix(T, N);
	this->densities = CallocDoubleMatrix(N, T);
	this->proba = (double*) Calloc(N, double);
	this->gamma = CallocDoubleMatrix(N, T);
	this->states_prev = (bool*) Calloc(T, bool);
	this->sumgamma = (double*) Calloc(N, double);
	this->sumxi = CallocDoubleMatrix(N, N);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
	this->sumdiff_state_last = 0;
// 	this->num_nonzero_A_into_state = (int*) Calloc(N, int);
// 	this->index_nonzero_A_into_state = allocIntMatrix(N, N);
// 	this->transition_cutoff = 1e-10;
// 	this->sparsity_cutoff = 0.0;

}

ScaleHMM::ScaleHMM(int T, int N, int Nmod, int verbosity)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	//FILE_LOG(logDEBUG2) << "Initializing multivariate ScaleHMM";
	this->xvariate = MULTIVARIATE;
	this->verbosity = verbosity;
	this->T = T;
	this->N = N;
	this->A = CallocDoubleMatrix(N, N);
	this->scalefactoralpha = (double*) Calloc(T, double);
	this->scalealpha = CallocDoubleMatrix(T, N);
	this->scalebeta = CallocDoubleMatrix(T, N);
	this->densities = CallocDoubleMatrix(N, T);
// 	this->tdensities = CallocDoubleMatrix(T, N);
	this->proba = (double*) Calloc(N, double);
	this->gamma = CallocDoubleMatrix(N, T);
	this->sumgamma = (double*) Calloc(N, double);
	this->sumxi = CallocDoubleMatrix(N, N);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
// 	this->sumdiff_state_last = 0;
	this->Nmod = Nmod;
// 	this->num_nonzero_A_into_state = (int*) Calloc(N, int);
// 	this->index_nonzero_A_into_state = allocIntMatrix(N, N);
// 	this->transition_cutoff = 1e-10;
// 	this->sparsity_cutoff = 0.7;

}

ScaleHMM::~ScaleHMM()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	FreeDoubleMatrix(this->A, this->N);
	Free(this->scalefactoralpha);
	FreeDoubleMatrix(this->scalealpha, this->T);
	FreeDoubleMatrix(this->scalebeta, this->T);
	FreeDoubleMatrix(this->densities, this->N);
// 	if (this->xvariate==MULTIVARIATE)
// 	{
// 		FreeDoubleMatrix(this->tdensities, this->T);
// 	}
	FreeDoubleMatrix(this->gamma, this->N);
	FreeDoubleMatrix(this->sumxi, this->N);
	Free(this->proba);
	Free(this->sumgamma);

	for (int iN=0; iN<this->N; iN++)
	{
		//FILE_LOG(logDEBUG1) << "Deleting density functions"; 
		delete this->densityFunctions[iN];
	}
// 	Free(this->num_nonzero_A_into_state);
// 	FreeIntMatrix(this->index_nonzero_A_into_state, this->N);
}

// Methods ----------------------------------------------------
void ScaleHMM::initialize_transition_probs(double* initial_A, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			for (int jN=0; jN<this->N; jN++)
			{
				// convert from vector to matrix representation
				this->A[jN][iN] = initial_A[iN*this->N + jN];
			}
		}
	}
	else
	{
		double self = 0.9;
// 		self = 1.0 / this->N; // set to uniform
		double other = (1.0 - self) / (this->N - 1.0);
		for (int iN=0; iN<this->N; iN++)
		{
			for (int jN=0; jN<this->N; jN++)
			{
				if (iN == jN)
					this->A[iN][jN] = self;
				else
					this->A[iN][jN] = other;
				// Save value to initial A
				initial_A[jN*this->N + iN] = this->A[iN][jN];
			}
		}
	}
	
// 	// Initialize sparsity exploit such that no sparsity exploit is done in first iteration
// 	for (int iN=0; iN<this->N; iN++)
// 	{
// 		this->num_nonzero_A_into_state[iN] = this->N;
// 	}
	

}

void ScaleHMM::initialize_proba(double* initial_proba, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = initial_proba[iN];
		}
	}
	else
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = (double)1/this->N;
			// Save value to initial proba
			initial_proba[iN] = this->proba[iN];
		}
	}

}

void ScaleHMM::baumWelch(int* maxiter, int* maxtime, double* eps)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;

	double logPold = -INFINITY;
	double logPnew;

	// Parallelization settings
// 	#ifdef _OPENMP
// 	omp_set_nested(1);
// 	#endif
	
	// measuring the time
	this->baumWelchStartTime_sec = time(NULL);

	if (this->xvariate == UNIVARIATE)
	{
		//FILE_LOG(logINFO) << "";
// 		Rprintf("\n");
		//FILE_LOG(logINFO) << "INITIAL PARAMETERS";
// 		Rprintf("INITIAL PARAMETERS\n");
		this->print_uni_params();
		this->print_uni_iteration(0);
	}
	else if (this->xvariate == MULTIVARIATE)
	{
		this->print_multi_iteration(0);
		//FILE_LOG(logDEBUG2) << "Calling calc_densities() from baumWelch()";
		//FILE_LOG(logINFO) << "Precomputing densities ...";
		if (this->verbosity>=1) Rprintf("HMM: Precomputing densities ...\n");
		try { this->calc_densities(); } catch(...) { throw; }
		this->print_multi_iteration(0);
		// Print densities
// 		int bs = 100;
// 		char buffer [100];
// 		int cx;
// 		for (int t=0; t<this->T; t++)
// 		{
// 			cx = 0;
// 			cx += snprintf(buffer+cx, bs-cx, "t=%d\t", t);
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				cx += snprintf(buffer+cx, bs-cx, "%.6f\t", this->densities[iN][t]);
// 			}
// 			//FILE_LOG(logDEBUG) << buffer;
// 		}
				
	}

	R_CheckUserInterrupt();
	// Do the Baum-Welch
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	int iteration = 0;
	while (((this->baumWelchTime_real < *maxtime) or (*maxtime < 0)) and ((iteration < *maxiter) or (*maxiter < 0)))
	{

		iteration++;
		
		if (this->xvariate == UNIVARIATE)
		{
			//FILE_LOG(logDEBUG1) << "Calling calc_densities() from baumWelch()";
			try { this->calc_densities(); } catch(...) { throw; }
			R_CheckUserInterrupt();
		}

		//FILE_LOG(logDEBUG1) << "Calling forward() from baumWelch()";
		try { this->forward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling backward() from baumWelch()";
		try { this->backward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_loglikelihood() from baumWelch()";
		this->calc_loglikelihood();
		logPnew = this->logP;
		if(std::isnan(logPnew))
		{
			//FILE_LOG(logERROR) << "logPnew = " << logPnew;
			throw nan_detected;
		}
		this->dlogP = logPnew - logPold;

		//FILE_LOG(logDEBUG1) << "Calling calc_sumxi() from baumWelch()";
		this->calc_sumxi();
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_sumgamma() from baumWelch()";
		this->calc_sumgamma();
		R_CheckUserInterrupt();

		if (this->xvariate == UNIVARIATE)
		{
// 			clock_t clocktime = clock(), dtime;
			// difference in state assignments
			//FILE_LOG(logDEBUG1) << "Calculating differences in state assignments in baumWelch()";
			int state_last = 0;
			int state_last_old = 0;
			int statesum = 0;
			for (int t=0; t<this->T; t++)
			{
				state_last_old = this->states_prev[t];
				if (this->gamma[this->N-1][t]>0.5)
				{
					state_last = 1;
				}
				this->states_prev[t] = state_last;
				statesum += std::abs(state_last-state_last_old);
				state_last = 0;
				state_last_old = 0;
			}
			this->sumdiff_state_last = statesum;

		}
		R_CheckUserInterrupt();

		// Print information about current iteration
		if (this->xvariate == UNIVARIATE)
		{
			this->print_uni_iteration(iteration);
		}
		else if (this->xvariate == MULTIVARIATE)
		{
			this->print_multi_iteration(iteration);
		}

		// Check convergence
		if(this->dlogP < *eps) //it has converged
		{
			//FILE_LOG(logINFO) << "Convergence reached!\n";
			if (this->verbosity>=1) Rprintf("HMM: Convergence reached!\n");
			if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
			break;
		} else {// not converged
			this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
			if (iteration == *maxiter)
			{
				//FILE_LOG(logINFO) << "Maximum number of iterations reached!";
				if (this->verbosity>=1) Rprintf("HMM: Maximum number of iterations reached!\n");
				if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
				break;
			}
			else if ((this->baumWelchTime_real >= *maxtime) and (*maxtime >= 0))
			{
				//FILE_LOG(logINFO) << "Exceeded maximum time!";
				if (this->verbosity>=1) Rprintf("HMM: Exceeded maximum time!\n");
				if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
				break;
			}
			logPold = logPnew;
		}

// 		// Re-Initialize the sparsity exploit variables
// 		int nonzero_counter [this->N];
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			this->num_nonzero_A_into_state[iN] = 0;
// 			nonzero_counter[iN] = 0;
// 			for (int jN=0; jN<this->N; jN++)
// 			{
// 				this->index_nonzero_A_into_state[iN][jN] = 0;
// 			}
// 		}
		
		// Updating initial probabilities proba and transition matrix A
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = this->gamma[iN][0];
			//FILE_LOG(logDEBUG4) << "sumgamma["<<iN<<"] = " << sumgamma[iN];
			if (this->sumgamma[iN] == 0)
			{
				//FILE_LOG(logINFO) << "Not reestimating A["<<iN<<"][x] because sumgamma["<<iN<<"] = 0";
// 				Rprintf("Not reestimating A[%d][x] because sumgamma[%d] = 0\n", iN, iN);
			}
			else
			{
				for (int jN=0; jN<this->N; jN++)
				{
					//FILE_LOG(logDEBUG4) << "sumxi["<<iN<<"]["<<jN<<"] = " << sumxi[iN][jN];
					this->A[iN][jN] = this->sumxi[iN][jN] / this->sumgamma[iN];
// 					// This could give performance increase, but risks numerical instabilities
// 					if (this->logA[iN][jN] < log(1.0/(double)(this->T*10)))
// 					{
// 						this->logA[iN][jN] = -INFINITY;
// 						this->A[iN][jN] = 0;
// 					}
// 					// Save the indices of non-zero transitions for sparsity exploit. We also get numerical instabilities here.
// 					if (this->A[iN][jN] > this->transition_cutoff)
// 					{
// 						this->num_nonzero_A_into_state[jN]++;
// 						this->index_nonzero_A_into_state[jN][nonzero_counter[jN]] = iN;
// 						nonzero_counter[jN]++;
// 					}
					if (std::isnan(this->A[iN][jN]))
					{
						//FILE_LOG(logERROR) << "updating transition probabilities";
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
						//FILE_LOG(logERROR) << "sumxi["<<iN<<"]["<<jN<<"] = " << sumxi[iN][jN];
						//FILE_LOG(logERROR) << "sumgamma["<<iN<<"] = " << sumgamma[iN];
						throw nan_detected;
					}
				}
			}
		}

// 		// Check if sparsity exploit will be used in the next iteration
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			if (this->num_nonzero_A_into_state[iN] < this->sparsity_cutoff*this->N)
// 			{
// 				//FILE_LOG(logINFO) << "Will use sparsity exploit for state " << iN;
// 				if (this->num_nonzero_A_into_state[iN] == 0)
// 				{
// 					//FILE_LOG(logINFO) << "No non-zero elements into state " << iN;
// 				}
// 			}
// 		}

		if (this->xvariate == UNIVARIATE)
		{
			// Update the parameters of the distribution
// 			clock_t clocktime = clock(), dtime;
			#pragma omp parallel for
			for (int iN=0; iN<this->N; iN++)
			{
				this->densityFunctions[iN]->update(this->gamma[iN]);
			}
// 			dtime = clock() - clocktime;
// 			//FILE_LOG(logDEBUG) << "updating distributions: " << dtime << " clicks";
			R_CheckUserInterrupt();
		}

	} /* main loop end */
    
    
	//Print the last results
	if (this->xvariate == UNIVARIATE)
	{
		//FILE_LOG(logINFO) << "";
// 		Rprintf("\n");
		//FILE_LOG(logINFO) << "FINAL ESTIMATION RESULTS";
// 		Rprintf("FINAL ESTIMATION RESULTS\n");
		this->print_uni_params();
	}

	// Return values
	*maxiter = iteration;
	*eps = this->dlogP;
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	*maxtime = this->baumWelchTime_real;
}

void ScaleHMM::check_for_state_swap()
{
	// only check for state swap if we have 3 states
	if (this->N == 3)
	{

		std::vector<double> weights(this->N);
		std::vector<double> maxdens(this->N);
		std::vector<double> cutoff_logdens(this->N);
		std::vector<double> logdens_at_0(this->N);

		// calculate weights, maxdens and logdens at cutoff
		weights = this->calc_weights();
		maxdens[0] = weights[0];
		maxdens[1] = weights[1] * Max(this->densities[1], this->T);
		maxdens[2] = weights[2] * Max(this->densities[2], this->T);
		cutoff_logdens[0] = log(weights[0]) + this->densityFunctions[0]->getLogDensityAt(this->cutoff);
		cutoff_logdens[1] = log(weights[1]) + this->densityFunctions[1]->getLogDensityAt(this->cutoff);
		cutoff_logdens[2] = log(weights[2]) + this->densityFunctions[2]->getLogDensityAt(this->cutoff);
		logdens_at_0[0] = log(weights[0]) + this->densityFunctions[0]->getLogDensityAt(0);
		logdens_at_0[1] = log(weights[1]) + this->densityFunctions[1]->getLogDensityAt(0);
		logdens_at_0[2] = log(weights[2]) + this->densityFunctions[2]->getLogDensityAt(0);

		// get highest value where both distributions are non-zero
		double logdens_1, logdens_2;
		bool doswap = false;
		for (int i1=this->cutoff; i1>=0; i1--)
		{
			logdens_1 = log(weights[1]) + this->densityFunctions[1]->getLogDensityAt(i1);
			logdens_2 = log(weights[2]) + this->densityFunctions[2]->getLogDensityAt(i1);
// 			Rprintf("i1 = %d\n",i1);
// 			Rprintf("logdens_1 = %g\n",logdens_1);
// 			Rprintf("logdens_2 = %g\n",logdens_2);
			if (logdens_1 != logdens_2)
			{
				if (logdens_1 > logdens_2)
				{
					doswap = true;
				}
				break;
			}
		}

		//FILE_LOG(logINFO) << "mean(0) = "<<this->densityFunctions[0]->get_mean() << ", mean(1) = "<<this->densityFunctions[1]->get_mean() << ", mean(2) = "<<this->densityFunctions[2]->get_mean();
// 		Rprintf("mean(0) = %g, mean(1) = %g, mean(2) = %g\n", this->densityFunctions[0]->get_mean(), this->densityFunctions[1]->get_mean(), this->densityFunctions[2]->get_mean());
		//FILE_LOG(logINFO) << "weight(0) = "<<weights[0] << ", weight(1) = "<<weights[1] << ", weight(2) = "<<weights[2];
// 		Rprintf("weight(0) = %g, weight(1) = %g, weight(2) = %g\n", weights[0], weights[1], weights[2]);
		//FILE_LOG(logINFO) << "maxdens(0) = "<<maxdens[0] << ", maxdens(1) = "<<maxdens[1] << ", maxdens(2) = "<<maxdens[2];
// 		Rprintf("maxdens(0) = %g, maxdens(1) = %g, maxdens(2) = %g\n", maxdens[0], maxdens[1], maxdens[2]);
		//FILE_LOG(logINFO) << "logdensity at x = "<<this->cutoff <<": logdens(0) = "<<cutoff_logdens[0] << ", logdens(1) = "<<cutoff_logdens[1] << ", logdens(2) = "<<cutoff_logdens[2];
// 		Rprintf("logdensity at x = %d: logdens(0) = %g, logdens(1) = %g, logdens(2) = %g\n", this->cutoff, cutoff_logdens[0], cutoff_logdens[1], cutoff_logdens[2]);
		// Different methods for state swapping detection
		// 1) Compare means. Does not work for all datasets.
	// 	if (this->densityFunctions[1]->get_mean() > this->densityFunctions[2]->get_mean()) //states 1 and 2 need to be exchanged
		// 2) Compare density values at upper cutoff. Works for most datasets.
// 		if (cutoff_logdens[1] > cutoff_logdens[2]) //states 1 and 2 need to be exchanged
		// 3) Compare max(density values). Does not work for all datasets.
	// 	if (maxdens[1] < maxdens[2])
		// 4) Compare density values at 0. Does not work for all datasets.
// 		if (logdens_at_0[1] < logdens_at_0[2])
		// 5) Compare logdens at highest value where both distributions are non-zero.
		if (doswap)
		{
			//FILE_LOG(logINFO) << "...swapping states";
// 			Rprintf("...swapping states\n");
			NegativeBinomial *tempDens = new NegativeBinomial();
			tempDens->copy(this->densityFunctions[2]); // tempDens is densifunc[2]
			this->densityFunctions[2]->copy(this->densityFunctions[1]); 
			this->densityFunctions[1]->copy(tempDens); 
			delete tempDens;
			// swap proba
			double temp;
			temp=this->proba[1];
			this->proba[1]=this->proba[2];
			this->proba[2]=temp;
			// swap transition matrix
			temp = this->A[0][2];
			this->A[0][2] = this->A[0][1];
			A[0][1] = temp;
			temp = this->A[1][0];
			this->A[1][0] = this->A[2][0];
			A[2][0] = temp;
			temp = this->A[1][1];
			this->A[1][1] = this->A[2][2];
			A[2][2] = temp;
			temp = this->A[1][2];
			this->A[1][2] = this->A[2][1];
			A[2][1] = temp;
			// swap dens
			double * tempp;
			tempp = this->densities[1];
			this->densities[1] = this->densities[2];
			this->densities[2] = tempp;
			// recalculate other baum-welch stuff
			//FILE_LOG(logDEBUG1) << "Calling forward() from check_for_state_swap()";
			try { this->forward(); } catch(...) { throw; }
			R_CheckUserInterrupt();
			//FILE_LOG(logDEBUG1) << "Calling backward() from check_for_state_swap()";
			try { this->backward(); } catch(...) { throw; }
			R_CheckUserInterrupt();
			//FILE_LOG(logDEBUG1) << "Calling calc_sumxi() from check_for_state_swap()";
			this->calc_sumxi();
			R_CheckUserInterrupt();
			//FILE_LOG(logDEBUG1) << "Calling calc_sumgamma() from check_for_state_swap()";
			this->calc_sumgamma();
			R_CheckUserInterrupt();
			
			// recalculate weight, maxdens and logdens at cutoff
			weights = this->calc_weights();
			maxdens[0] = weights[0];
			maxdens[1] = weights[1] * Max(this->densities[1], this->T);
			maxdens[2] = weights[2] * Max(this->densities[2], this->T);
			cutoff_logdens[0] = log(weights[0]) + this->densityFunctions[0]->getLogDensityAt(this->cutoff);
			cutoff_logdens[1] = log(weights[1]) + this->densityFunctions[1]->getLogDensityAt(this->cutoff);
			cutoff_logdens[2] = log(weights[2]) + this->densityFunctions[2]->getLogDensityAt(this->cutoff);
			logdens_at_0[0] = log(weights[0]) + this->densityFunctions[0]->getLogDensityAt(0);
			logdens_at_0[1] = log(weights[1]) + this->densityFunctions[1]->getLogDensityAt(0);
			logdens_at_0[2] = log(weights[2]) + this->densityFunctions[2]->getLogDensityAt(0);

			//FILE_LOG(logINFO) << "mean(0) = "<<this->densityFunctions[0]->get_mean() << ", mean(1) = "<<this->densityFunctions[1]->get_mean() << ", mean(2) = "<<this->densityFunctions[2]->get_mean();
// 			Rprintf("mean(0) = %g, mean(1) = %g, mean(2) = %g\n", this->densityFunctions[0]->get_mean(), this->densityFunctions[1]->get_mean(), this->densityFunctions[2]->get_mean());
			//FILE_LOG(logINFO) << "weight(0) = "<<weights[0] << ", weight(1) = "<<weights[1] << ", weight(2) = "<<weights[2];
// 			Rprintf("weight(0) = %g, weight(1) = %g, weight(2) = %g\n", weights[0], weights[1], weights[2]);
			//FILE_LOG(logINFO) << "maxdens(0) = "<<maxdens[0] << ", maxdens(1) = "<<maxdens[1] << ", maxdens(2) = "<<maxdens[2];
// 			Rprintf("maxdens(0) = %g, maxdens(1) = %g, maxdens(2) = %g\n", maxdens[0], maxdens[1], maxdens[2]);
			//FILE_LOG(logINFO) << "logdensity at x = "<<this->cutoff <<": logdens(0) = "<<cutoff_logdens[0] << ", logdens(1) = "<<cutoff_logdens[1] << ", logdens(2) = "<<cutoff_logdens[2];
// 			Rprintf("logdensity at x = %d: logdens(0) = %g, logdens(1) = %g, logdens(2) = %g\n", this->cutoff, cutoff_logdens[0], cutoff_logdens[1], cutoff_logdens[2]);
		}
	}
}

std::vector<double> ScaleHMM::calc_weights()
{
	std::vector<double> weights(this->N);
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->T; t++)
		{
			sum_over_gammas_per_state += this->gamma[iN][t];
		}
		weights[iN] = sum_over_gammas_per_state / this->T;
	}
	return(weights);
}

void ScaleHMM::calc_weights(double* weights)
{
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->T; t++)
		{
			sum_over_gammas_per_state += this->gamma[iN][t];
		}
		weights[iN] = sum_over_gammas_per_state / this->T;
	}
}

// Getters and Setters ----------------------------------------
int ScaleHMM::get_N()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->N);
}

int ScaleHMM::get_T()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->T);
}

void ScaleHMM::get_posteriors(double** post)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int iN=0; iN<this->N; iN++)
	{
		for (int t=0; t<this->T; t++)
		{
			post[iN][t] = this->gamma[iN][t];
		}
	}
}

double ScaleHMM::get_posterior(int iN, int t)
{
	//FILE_LOG(logDEBUG4) << __PRETTY_FUNCTION__;
	return(this->gamma[iN][t]);
}

double ScaleHMM::get_density(int iN, int t)
{
	//FILE_LOG(logDEBUG3) << __PRETTY_FUNCTION__;
	return(this->densities[iN][t]);
}

double ScaleHMM::get_proba(int i)
{
	return( this->proba[i] );
}

double ScaleHMM::get_A(int i, int j)
{
	return( this->A[i][j] );
}

double ScaleHMM::get_logP()
{
	return( this->logP );
}

void ScaleHMM::set_cutoff(int cutoff)
{
	this->cutoff = cutoff;
}

// Private ====================================================
// Methods ----------------------------------------------------
void ScaleHMM::forward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

		std::vector<double> alpha(this->N);
		// Initialization
		this->scalefactoralpha[0] = 0.0;
		for (int iN=0; iN<this->N; iN++)
		{
			alpha[iN] = this->proba[iN] * this->densities[iN][0];
			//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
			this->scalefactoralpha[0] += alpha[iN];
		}
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<0<<"] = " << scalefactoralpha[0];
		for (int iN=0; iN<this->N; iN++)
		{
			this->scalealpha[0][iN] = alpha[iN] / this->scalefactoralpha[0];
			//FILE_LOG(logDEBUG4) << "scalealpha["<<0<<"]["<<iN<<"] = " << scalealpha[0][iN];
		}
		// Induction
		for (int t=1; t<this->T; t++)
		{
			this->scalefactoralpha[t] = 0.0;
			for (int iN=0; iN<this->N; iN++)
			{
				double helpsum = 0.0;
				for (int jN=0; jN<this->N; jN++)
				{
					helpsum += this->scalealpha[t-1][jN] * this->A[jN][iN];
				}
				alpha[iN] = helpsum * this->densities[iN][t];
				//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
				this->scalefactoralpha[t] += alpha[iN];
			}
			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalealpha[t][iN] = alpha[iN] / this->scalefactoralpha[t];
				//FILE_LOG(logDEBUG4) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
				if(std::isnan(this->scalealpha[t][iN]))
				{
					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					for (int jN=0; jN<this->N; jN++)
					{
						//FILE_LOG(logERROR) << "scalealpha["<<t-1<<"]["<<jN<<"] = " << scalealpha[t-1][jN];
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
					}
					//FILE_LOG(logERROR) << "scalefactoralpha["<<t<<"] = "<<scalefactoralpha[t] << ", densities = "<<densities[iN][t];
					//FILE_LOG(logERROR) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		std::vector<double> alpha(this->N);
// 		// Initialization
// 		this->scalefactoralpha[0] = 0.0;
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			alpha[iN] = this->proba[iN] * this->tdensities[0][iN];
// 			//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
// 			this->scalefactoralpha[0] += alpha[iN];
// 		}
// 		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<0<<"] = " << scalefactoralpha[0];
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			this->scalealpha[0][iN] = alpha[iN] / this->scalefactoralpha[0];
// 			//FILE_LOG(logDEBUG4) << "scalealpha["<<0<<"]["<<iN<<"] = " << scalealpha[0][iN];
// 		}
// 		// Induction
// 		for (int t=1; t<this->T; t++)
// 		{
// 			this->scalefactoralpha[t] = 0.0;
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				double helpsum = 0.0;
// 				for (int jN=0; jN<this->N; jN++)
// 				{
// 					helpsum += this->scalealpha[t-1][jN] * this->A[jN][iN];
// 				}
// 				alpha[iN] = helpsum * this->tdensities[t][iN];
// 				//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
// 				this->scalefactoralpha[t] += alpha[iN];
// 			}
// 			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->scalealpha[t][iN] = alpha[iN] / this->scalefactoralpha[t];
// 				//FILE_LOG(logDEBUG4) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
// 				if(std::isnan(this->scalealpha[t][iN]))
// 				{
// 					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
// 					for (int jN=0; jN<this->N; jN++)
// 					{
// 						//FILE_LOG(logERROR) << "scalealpha["<<t-1<<"]["<<jN<<"] = " << scalealpha[t-1][jN];
// 						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
// 					}
// 					//FILE_LOG(logERROR) << "scalefactoralpha["<<t<<"] = "<<scalefactoralpha[t] << ", tdensities = "<<tdensities[t][iN];
// 					//FILE_LOG(logERROR) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
// 					throw nan_detected;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "forward(): " << dtime << " clicks";
}

void ScaleHMM::backward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

		std::vector<double> beta(this->N);
		// Initialization
		for (int iN=0; iN<this->N; iN++)
		{
			beta[iN] = 1.0;
			//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
		}
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<T-1<<"] = " << scalefactoralpha[T-1];
		for (int iN=0; iN<this->N; iN++)
		{
			this->scalebeta[T-1][iN] = beta[iN] / this->scalefactoralpha[T-1];
			//FILE_LOG(logDEBUG4) << "scalebeta["<<T-1<<"]["<<iN<<"] = " << scalebeta[T-1][iN];
		}
		// Induction
		for (int t=this->T-2; t>=0; t--)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				beta[iN] = 0.0;
				for(int jN=0; jN<this->N; jN++)
				{
					beta[iN] += this->A[iN][jN] * this->densities[jN][t+1] * this->scalebeta[t+1][jN];
				}
				//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
			}
			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalebeta[t][iN] = beta[iN] / this->scalefactoralpha[t];
				//FILE_LOG(logDEBUG4) << "scalebeta["<<t<<"]["<<iN<<"] = " << scalebeta[t][iN];
				if (std::isnan(this->scalebeta[t][iN]))
				{
					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					for (int jN=0; jN<this->N; jN++)
					{
						//FILE_LOG(logERROR) << "scalebeta["<<jN<<"]["<<t+1<<"] = " << scalebeta[t+1][jN];
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
						//FILE_LOG(logERROR) << "densities["<<jN<<"]["<<t+1<<"] = " << densities[jN][t+1];
					}
					//FILE_LOG(logERROR) << "this->scalefactoralpha[t]["<<t<<"] = "<<this->scalefactoralpha[t] << ", densities = "<<densities[iN][t];
					//FILE_LOG(logERROR) << "scalebeta["<<iN<<"]["<<t<<"] = " << scalebeta[t][iN];
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		std::vector<double> beta(this->N);
// 		// Initialization
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			beta[iN] = 1.0;
// 			//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
// 		}
// 		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<T-1<<"] = " << scalefactoralpha[T-1];
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			this->scalebeta[T-1][iN] = beta[iN] / this->scalefactoralpha[T-1];
// 			//FILE_LOG(logDEBUG4) << "scalebeta["<<T-1<<"]["<<iN<<"] = " << scalebeta[T-1][iN];
// 		}
// 		// Induction
// 		for (int t=this->T-2; t>=0; t--)
// 		{
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				//FILE_LOG(logDEBUG4) << "Calculating backward variable for state " << iN;
// 				beta[iN] = 0.0;
// 				for(int jN=0; jN<this->N; jN++)
// 				{
// 					beta[iN] += this->A[iN][jN] * this->tdensities[t+1][jN] * this->scalebeta[t+1][jN];
// 				}
// 			}
// 			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->scalebeta[t][iN] = beta[iN] / this->scalefactoralpha[t];
// 				//FILE_LOG(logDEBUG4) << "scalebeta["<<t<<"]["<<iN<<"] = " << scalebeta[t][iN];
// 				if (std::isnan(this->scalebeta[t][iN]))
// 				{
// 					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
// 					for (int jN=0; jN<this->N; jN++)
// 					{
// 						//FILE_LOG(logERROR) << "scalebeta["<<jN<<"]["<<t+1<<"] = " << scalebeta[t+1][jN];
// 						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
// 						//FILE_LOG(logERROR) << "tdensities["<<t+1<<"]["<<jN<<"] = " << tdensities[t+1][jN];
// 					}
// 					//FILE_LOG(logERROR) << "this->scalefactoralpha[t]["<<t<<"] = "<<this->scalefactoralpha[t] << ", tdensities = "<<densities[t][iN];
// 					//FILE_LOG(logERROR) << "scalebeta["<<iN<<"]["<<t<<"] = " << scalebeta[t][iN];
// 					throw nan_detected;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "backward(): " << dtime << " clicks";
}

void ScaleHMM::calc_sumgamma()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	// Initialize the sumgamma
	for (int iN=0; iN<this->N; iN++)
	{
		this->sumgamma[iN] = 0.0;
	}

	// Compute the gammas (posteriors) and sumgamma
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		for (int t=0; t<this->T; t++)
		{
			this->gamma[iN][t] = this->scalealpha[t][iN] * this->scalebeta[t][iN] * this->scalefactoralpha[t];
			this->sumgamma[iN] += this->gamma[iN][t];
		}
	}
	// Subtract the last value because sumgamma goes only until T-1 and we computed until T to get also loggamma at T
	for (int iN=0; iN<this->N; iN++)
	{
		this->sumgamma[iN] -= this->gamma[iN][T-1];
	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_sumgamma(): " << dtime << " clicks";
}

void ScaleHMM::calc_sumxi()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	double xi;
	// Initialize the sumxi
	for (int iN=0; iN<this->N; iN++)
	{
		for (int jN=0; jN<this->N; jN++)
		{
			this->sumxi[iN][jN] = 0.0;
		}
	}	

// 	if (this->xvariate==UNIVARIATE)
// 	{

		#pragma omp parallel for
		for (int iN=0; iN<this->N; iN++)
		{
			//FILE_LOG(logDEBUG3) << "Calculating sumxi["<<iN<<"][jN]";
			for (int t=0; t<this->T-1; t++)
			{
				for (int jN=0; jN<this->N; jN++)
				{
					xi = this->scalealpha[t][iN] * this->A[iN][jN] * this->densities[jN][t+1] * this->scalebeta[t+1][jN];
					this->sumxi[iN][jN] += xi;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		#pragma omp parallel for
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			//FILE_LOG(logDEBUG3) << "Calculating sumxi["<<iN<<"][jN]";
// 			for (int t=0; t<this->T-1; t++)
// 			{
// 				for (int jN=0; jN<this->N; jN++)
// 				{
// 					xi = this->scalealpha[t][iN] * this->A[iN][jN] * this->tdensities[t+1][jN] * this->scalebeta[t+1][jN];
// 					this->sumxi[iN][jN] += xi;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_sumxi(): " << dtime << " clicks";
}

void ScaleHMM::calc_loglikelihood()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	this->logP = 0.0;
	for (int t=0; t<this->T; t++)
	{
		this->logP += log(this->scalefactoralpha[t]);
	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_loglikelihood(): " << dtime << " clicks";
}

void ScaleHMM::calc_densities()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;
	// Errors thrown inside a #pragma must be handled inside the thread
	std::vector<bool> nan_encountered(this->N);
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		//FILE_LOG(logDEBUG3) << "Calculating densities for state " << iN;
		try
		{
			this->densityFunctions[iN]->calc_densities(this->densities[iN]);
		}
		catch(std::exception& e)
		{
			if (strcmp(e.what(),"nan detected")==0) { nan_encountered[iN]=true; }
			else { throw; }
		}
	}
	for (int iN=0; iN<this->N; iN++)
	{
		if (nan_encountered[iN]==true)
		{
			throw nan_detected;
		}
	}

	// Check if the density for all states is close to zero and correct to prevent NaNs
// 	double zero_cutoff = 1.18e-37; // 32-bit precision is 1.18e-38
	double zero_cutoff = 2.23e-307; // 64-bit precision is 2.23e-308
	std::vector<double> temp(this->N);
	// t=0
	for (int iN=0; iN<this->N; iN++)
	{
		temp[iN] = this->densities[iN][0];
	}
	if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->densities[iN][0] = zero_cutoff;
		}
	}
	// t>0
	for (int t=1; t<this->T; t++)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			temp[iN] = this->densities[iN][t];
		}
		if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				this->densities[iN][t] = this->densities[iN][t-1];
			}
		}
	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_densities(): " << dtime << " clicks";
}

void ScaleHMM::print_uni_iteration(int iteration)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->verbosity>=1)
	{
		this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
		int bs = 106;
		char buffer [106];
		if (iteration % 20 == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%19s%d%15s", "Iteration", "log(P)", "dlog(P)", "Diff in state ",this->N-1, "Time in sec");
			//FILE_LOG(logITERATION) << buffer;
			Rprintf("%s\n", buffer);
		}
		if (iteration == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%20s%*d", "0", "-inf", "-", "-", 15, this->baumWelchTime_real);
		}
		else if (iteration == 1)
		{
			snprintf(buffer, bs, "%*d%*f%20s%*d%*d", 10, iteration, 20, this->logP, "inf", 20, this->sumdiff_state_last, 15, this->baumWelchTime_real);
		}
		else
		{
			snprintf(buffer, bs, "%*d%*f%*f%*d%*d", 10, iteration, 20, this->logP, 20, this->dlogP, 20, this->sumdiff_state_last, 15, this->baumWelchTime_real);
		}
		//FILE_LOG(logITERATION) << buffer;
		Rprintf("%s\n", buffer);

		// Flush Rprintf statements to R console
		R_FlushConsole();
	}
}

void ScaleHMM::print_multi_iteration(int iteration)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->verbosity>=1)
	{
		this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
		int bs = 86;
		char buffer [86];
		if (iteration % 20 == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%15s", "Iteration", "log(P)", "dlog(P)", "Time in sec");
			//FILE_LOG(logITERATION) << buffer;
			Rprintf("%s\n", buffer);
		}
		if (iteration == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%*d", "0", "-inf", "-", 15, this->baumWelchTime_real);
		}
		else if (iteration == 1)
		{
			snprintf(buffer, bs, "%*d%*f%20s%*d", 10, iteration, 20, this->logP, "inf", 15, this->baumWelchTime_real);
		}
		else
		{
			snprintf(buffer, bs, "%*d%*f%*f%*d", 10, iteration, 20, this->logP, 20, this->dlogP, 15, this->baumWelchTime_real);
		}
		//FILE_LOG(logITERATION) << buffer;
		Rprintf("%s\n", buffer);

		// Flush Rprintf statements to R console
		R_FlushConsole();
	}
}

void ScaleHMM::print_uni_params()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->verbosity>=2)
	{
		int bs = 82;
		char buffer [82];
		int cx;
		snprintf(buffer, bs, " -------------------------------------------------------------------------------");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		snprintf(buffer, bs, "|%80s", "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		// print loglik
		snprintf(buffer, bs, "| log(P) = %*.6f%54s", 16, this->logP, "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		snprintf(buffer, bs, "|%80s", "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		// print initial probabilities
		cx = snprintf(buffer, bs, "|%7s", "");
		for (int iN=0; iN<this->N; iN++)
		{
			cx += snprintf(buffer+cx, bs-cx, "proba[%d] = %.6f    ", iN, this->proba[iN]);
		}
		cx += snprintf(buffer+cx, bs-cx, "   |");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		snprintf(buffer, bs, "|%80s", "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		// print transition probabilities
		for (int iN=0; iN<this->N; iN++)
		{
			cx = snprintf(buffer, bs, "|%7s", "");
			for (int jN=0; jN<this->N; jN++)
			{
				cx += snprintf(buffer+cx, bs-cx, "A[%d][%d] = %.6f    ", iN, jN, this->A[iN][jN]);
			}
			cx += snprintf(buffer+cx, bs-cx, "      |");
			//FILE_LOG(logINFO) << buffer;
	 		Rprintf("%s\n", buffer);
		}
		// print emission parameters
		snprintf(buffer, bs, "|%80s", "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		for (int iN=0; iN<this->N; iN++)
		{
			if (iN == 1)
			{
				snprintf(buffer, bs, "| unmodified component%59s", "|");
				//FILE_LOG(logINFO) << buffer;
	 			Rprintf("%s\n", buffer);
			}
			if (iN == 2)
			{
				snprintf(buffer, bs, "| modified component%61s", "|");
				//FILE_LOG(logINFO) << buffer;
	 			Rprintf("%s\n", buffer);
			}
			if (this->densityFunctions[iN]->get_name() == NEGATIVE_BINOMIAL)
			{
				NegativeBinomial* temp = (NegativeBinomial*)this->densityFunctions[iN];
				double curR = temp->get_size();
				double curP = temp->get_prob();
				double curMean = temp->get_mean();
				double curVar = temp->get_variance();
				snprintf(buffer, bs, "| r = %*.6f, p = %*.6f, mean = %*.2f, var = %*.2f%20s", 9, curR, 9, curP, 6, curMean, 8, curVar, "|");
				//FILE_LOG(logINFO) << buffer;
	 			Rprintf("%s\n", buffer);
			}
		}
		
		snprintf(buffer, bs, "|%80s", "|");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		snprintf(buffer, bs, " -------------------------------------------------------------------------------");
		//FILE_LOG(logINFO) << buffer;
	 	Rprintf("%s\n", buffer);
		//FILE_LOG(logINFO) << "";
	 	Rprintf("\n");

		// Flush Rprintf statements to R console
	 	R_FlushConsole();
	}
}

