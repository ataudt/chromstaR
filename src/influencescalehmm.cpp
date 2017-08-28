
#include "influencescalehmm.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------


// InfluenceScaleHMM::InfluenceScaleHMM(int T, int N, int verbosity)
// {
// 	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	//FILE_LOG(logDEBUG2) << "Initializing univariate InfluenceScaleHMM";
// 	this->xvariate = UNIVARIATE;
// 	this->verbosity = verbosity;
// 	this->T = T;
// 	this->N = N;
// 	this->A = alloc4Ddouble(Nmod,Nmod,N, N);
// 	this->scalefactoralpha = CallocDoubleMatrix(Nmod,T);
// 	this->scalealpha = alloc3Ddouble(Nmod,T, N);
// 	this->scalebeta = alloc3Ddouble(Nmod,T, N);
// 	this->densities = alloc3Ddouble(Nmod,N, T);
// 	this->proba = CallocDoubleMatrix(Nmod,N);
// 	this->gamma = alloc4Ddouble(Nmod,Nmod,N, T);
// 	this->states_prev = (bool*) Calloc(T, bool);
// 	this->sumgamma = alloc3Ddouble(Nmod,Nmod,N);
// 	//this->sumxi = CallocDoubleMatrix(N, N);
// 	this->logP = -INFINITY;
// 	this->dlogP = INFINITY;
// 	this->sumdiff_state_last = 0;
// // 	this->num_nonzero_A_into_state = (int*) Calloc(N, int);
// // 	this->index_nonzero_A_into_state = allocIntMatrix(N, N);
// // 	this->transition_cutoff = 1e-10;
// // 	this->sparsity_cutoff = 0.0;
//
// }


InfluenceScaleHMM::InfluenceScaleHMM(int T, int N, int Nmod, int verbosity)
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	//FILE_LOG(logDEBUG2) << "Initializing multivariate InfluenceScaleHMM";
	this->xvariate = MULTIVARIATE;
	this->verbosity = verbosity;
	this->T = T;
	this->N = N;
	this->A = alloc4Ddouble(Nmod,Nmod,N, N);
	this->scalefactoralpha = CallocDoubleMatrix(Nmod,T);
	this->scalealpha = alloc3Ddouble(Nmod,T, N);
	this->scalebeta = alloc3Ddouble(Nmod,T, N);
	this->densities = alloc3Ddouble(Nmod,N, T);
	//this->tdensities = CallocDoubleMatrix(T, N);
	this->proba = CallocDoubleMatrix(Nmod,N);
	this->gamma = alloc3Ddouble(Nmod,N, T);
	//this->sumgamma = CallocDoubleMatrix(Nmod,N);
	this->sumxi = alloc4Ddouble(Nmod, Nmod, N, N);
	this->tiesumxi= CallocDoubleMatrix(Nmod, Nmod);
	this->tie_onestatus_sumxi= alloc3Ddouble(Nmod, Nmod, N);
	//this->tiesumxitotal=CallocDoubleMatrix(Nmod,Nmod);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
  // 	this->sumdiff_state_last = 0;
	this->Nmod = Nmod;
// 	this->num_nonzero_A_into_state = (int*) Calloc(N, int);
// 	this->index_nonzero_A_into_state = allocIntMatrix(N, N);
// 	this->transition_cutoff = 1e-10;
// 	this->sparsity_cutoff = 0.7;
	this->influence= alloc4Ddouble(Nmod,Nmod,N,N);
	this->tiestrength= CallocDoubleMatrix(Nmod,Nmod);
	this->weights= CallocDoubleMatrix(Nmod,N);



}


InfluenceScaleHMM::~InfluenceScaleHMM()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	free4Ddouble(this->A,this->Nmod, this->Nmod, this->N);
	FreeDoubleMatrix(this->scalefactoralpha, this->Nmod);
	free3Ddouble(this->scalealpha,this->Nmod, this->T);
	free3Ddouble(this->scalebeta, this->Nmod,this->T);
	free3Ddouble(this->densities, this->Nmod,this->N);
// 	if (this->xvariate==MULTIVARIATE)
// 	{
// 		FreeDoubleMatrix(this->tdensities, this->T);
// 	}
	free3Ddouble(this->gamma,this->Nmod, this->Nmod);
	FreeDoubleMatrix(this->weights, this->Nmod);
	free4Ddouble(this->sumxi, this->Nmod, this->Nmod, this->N);
	free3Ddouble(this->tie_onestatus_sumxi, this->Nmod, this->Nmod);
	FreeDoubleMatrix(this->proba, this->Nmod);
	free4Ddouble(this->influence, this->Nmod, this->Nmod, this->N);
	FreeDoubleMatrix(this->scalefactoralpha, this->Nmod);

	for (int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++) {
		for (int iN=0; iN<this->N; iN++)
		{
			//FILE_LOG(logDEBUG1) << "Deleting density functions";
			delete this->densityFunctions[c1Nmod][iN];
		}
	}
	// for (int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
	// {
	// 	//FILE_LOG(logDEBUG1) << "Deleting density functions";
	// 	delete this->tiesumxitotal[c1Nmod];
	// }
// 	Free(this->num_nonzero_A_into_state);
// 	FreeIntMatrix(this->index_nonzero_A_into_state, this->N);
}


//methods

void InfluenceScaleHMM::initialize_transition_probs(double* initial_A, bool use_initial_params)
{

	if (use_initial_params)
	{
 		int i1=0 ;
		for (int jN=0; jN<this->N; jN++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
				{
					for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
					{
						// convert from vector to matrix representation
						this->A[c1Nmod][c2Nmod][iN][jN] = initial_A[i1];
						// Rprintf(" c1Nmod=%d c2Nmod=%d iN=%d jN=%d A=%g \n", c1Nmod, c2Nmod, iN, jN, this->A[c1Nmod][c2Nmod][iN][jN] );
						 i1++;

					}
				}
			}
		}

	}
	else
	{
		double self = 0.9;
		double other = (1.0 - self) / ((this->N - 1.0));

		int i1=0 ;
		for (int jN=0; jN<this->N; jN++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
				{
					for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
					{

						if (iN == jN )
							this->A[c1Nmod][c2Nmod][iN][jN] = self;
						else
							this->A[c1Nmod][c2Nmod][iN][jN] = other;
						// Save value to initial A
						initial_A[i1] = this->A[c1Nmod][c2Nmod][iN][jN];

						i1++;
					}
				}
			}
		}
	}

// 	// Initialize sparsity exploit such that no sparsity exploit is done in first iteration
// 	for (int iN=0; iN<this->N; iN++)
// 	{
// 		this->num_nonzero_A_into_state[iN] = this->N;
// 	}


}


void InfluenceScaleHMM::initialize_proba(double* initial_proba, bool use_initial_params)
{

	if (use_initial_params)
	{
		int i1=0;
		for (int iN=0; iN<this->N; iN++)
		{
			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
			{
				this->proba[c1Nmod][iN] = initial_proba[i1];
				i1++;
			}
		}
	}
	else
	{
		int i1=0;
		for (int iN=0; iN<this->N; iN++)
		{
			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
			{
				this->proba[c1Nmod][iN] = (double)1/((this->N));
				// Save value to initial proba
				initial_proba[i1] = this->proba[c1Nmod][iN];
				i1++;
			}
		}
	}

}




//--- LUISA these initialization functions are not based on an initialized vector like it is done for proba and transisions. Should I do it like that?
//-- As off now every value is zero.

void InfluenceScaleHMM::initialize_influence(){

	for(int c1Nmod=0; c1Nmod<this->Nmod ; c1Nmod++)
	{
		for(int c2Nmod=0; c2Nmod<this->Nmod ; c2Nmod++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				for (int jN=0; jN<this->N; jN++)
				{
					this->influence[c1Nmod][c2Nmod][iN][jN]= this->tiestrength[c1Nmod][c2Nmod] * this-> A[c1Nmod][c2Nmod][iN][jN];
				}
			}
		}
	}
}

//initial_tiestrength has dimension as initial_proba

void InfluenceScaleHMM::initialize_tiestrength(double* initial_tiestrength, bool use_initial_params){

	if (use_initial_params)
	{
		int i1=0;
		for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
		{
			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
			{
				this->tiestrength[c1Nmod][c2Nmod]= initial_tiestrength[i1];
				i1++;
			}
		}
	}
	else
	{
		int i1=0;
		for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
		{
			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
			{
				this->tiestrength[c1Nmod][c2Nmod]= double(1)/(this->Nmod);
				initial_tiestrength[i1] = this->tiestrength[c1Nmod][c2Nmod];
				i1++;
			}
		}
	}

}



void InfluenceScaleHMM::baumWelch(int* maxiter, int* maxtime, double* eps)
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}
	double logPold = -INFINITY;
	double logPnew;

	// Parallelization settings
// 	#ifdef _OPENMP
// 	omp_set_nested(1);
// 	#endif

	// measuring the time
	this->baumWelchStartTime_sec = time(NULL);

	//if (this->xvariate == UNIVARIATE)
	//{
		//FILE_LOG(logINFO) << "";
// 		Rprintf("\n");
		//FILE_LOG(logINFO) << "INITIAL PARAMETERS";
// 		Rprintf("INITIAL PARAMETERS\n");
	//	this->print_uni_params();
	//	this->print_uni_iteration(0);
	//}
	//else
 if (this->xvariate == MULTIVARIATE)
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

		this->calc_tiesumxi();
		R_CheckUserInterrupt();


//--LUISA UNIVARIATE ---should i let it like it is ???
//LUISA--comment
/**
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
		*/
		R_CheckUserInterrupt();

		// Print information about current iteration
		//if (this->xvariate == UNIVARIATE)
		//{
		//	this->print_uni_iteration(iteration);
		//}
		//else
		if (this->xvariate == MULTIVARIATE)
		{
			this->print_multi_iteration(iteration);
		}

		// Check convergence
		if(this->dlogP < *eps) //it has converged
		{
			//FILE_LOG(logINFO) << "Convergence reached!\n";
			if (this->verbosity>=1) Rprintf("HMM: Convergence reached!\n");
			//if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
			break;
		} else {// not converged
			this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
			if (iteration == *maxiter)
			{
				//FILE_LOG(logINFO) << "Maximum number of iterations reached!";
				if (this->verbosity>=1) Rprintf("HMM: Maximum number of iterations reached!\n");
				//if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
				break;
			}
			else if ((this->baumWelchTime_real >= *maxtime) and (*maxtime >= 0))
			{
				//FILE_LOG(logINFO) << "Exceeded maximum time!";
				if (this->verbosity>=1) Rprintf("HMM: Exceeded maximum time!\n");
				//if (this->xvariate == UNIVARIATE) this->check_for_state_swap();
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

		for(int c1Nmod=0; c1Nmod<this->Nmod ; c1Nmod++)
		{
			for(int c2Nmod=0; c2Nmod<this->Nmod ; c2Nmod++)
			{
				for (int iN=0; iN<this->N; iN++)
				{
					this->proba[c1Nmod][iN] = this->gamma[c1Nmod][iN][0];
					//FILE_LOG(logDEBUG4) << "sumgamma["<<iN<<"] = " << sumgamma[iN];
					//TODO : not sumgamma
					if (this->tie_onestatus_sumxi[c1Nmod][c2Nmod][iN] == 0)
					{
						//FILE_LOG(logINFO) << "Not reestimating A["<<iN<<"][x] because sumgamma["<<iN<<"] = 0";
		// 				Rprintf("Not reestimating A[%d][x] because sumgamma[%d] = 0\n", iN, iN);
					}
					else
					{
						for (int jN=0; jN<this->N; jN++)
						{
							//FILE_LOG(logDEBUG4) << "sumxi["<<iN<<"]["<<jN<<"] = " << sumxi[iN][jN];
							this->A[c1Nmod][c2Nmod][iN][jN] = this->sumxi[c1Nmod][c2Nmod][iN][jN] / this->tie_onestatus_sumxi[c1Nmod][c2Nmod][iN];

							//LUISA
							this->tiestrength[c1Nmod][c2Nmod]= this->tiesumxi[c1Nmod][c2Nmod]/this->tiesumxitotal[c2Nmod];

							this->influence[c1Nmod][c2Nmod][iN][jN]= this->tiestrength[c1Nmod][c2Nmod] * this-> A[c1Nmod][c2Nmod][iN][jN];

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
							if (std::isnan(this->A[c1Nmod][c2Nmod][iN][jN]))
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

//LUISA univariate so no touch for now
//comments
/*
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


	}
	*/
	/* main loop end */


	//Print the last results
	//if (this->xvariate == UNIVARIATE)
	//{
		//FILE_LOG(logINFO) << "";
// 		Rprintf("\n");
		//FILE_LOG(logINFO) << "FINAL ESTIMATION RESULTS";
// 		Rprintf("FINAL ESTIMATION RESULTS\n");
	//	this->print_uni_params();
//	}

	// Return values
	*maxiter = iteration;
	*eps = this->dlogP;
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	*maxtime = this->baumWelchTime_real;
}

//LUISA commented

/*
void InfluenceScaleHMM::check_for_state_swap()
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
//LUISA-----changing all to null for chains, or for all chains??
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
			temp = this->A[0][0][0][2];
			this->A[0][0][0][2] = this->A[0][0][0][1];
			A[0][1] = temp;
			temp = this->A[0][0][1][0];
			this->A[1][0] = this->A[0][0][2][0];
			A[0][0][2][0] = temp;
			temp = this->A[0][0][1][1];
			this->A[0][0][1][1] = this->A[0][0][2][2];
			A[0][0][2][2] = temp;
			temp = this->A[0][0][1][2];
			this->A[0][0][1][2] = this->A[0][0][2][1];
			A[0][0][2][1] = temp;
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

			this->calc_tiesumxi();
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

*/
}




std::vector< std::vector< double > > InfluenceScaleHMM::calc_weights()
{
//	std::vector<double> weights(this->N);
	std::vector< std::vector< double > > weights(this->N);
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++){
			for (int t=0; t<this->T; t++)
				{
					sum_over_gammas_per_state += this->gamma[c1Nmod][iN][t];
				}

				this->weights[c1Nmod][iN] = sum_over_gammas_per_state / this->T;
		}

	}
	return(weights);
}
//TODO see if needed
void InfluenceScaleHMM::calc_weights(double* weights)
{
	#pragma omp parallel for

	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for(int c1Nmod=0; c1Nmod<Nmod; c1Nmod++){
			for (int t=0; t<this->T; t++)
			{
				sum_over_gammas_per_state += this->gamma[c1Nmod][iN][t];
			}
			this->weights[c1Nmod][iN] = sum_over_gammas_per_state / this->T;
	}

	}
}

// Getters and Setters ----------------------------------------
int InfluenceScaleHMM::get_N()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	return(this->N);
}

int InfluenceScaleHMM::get_T()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	return(this->T);
}


//---LUISA--- change?????????--- problem
void InfluenceScaleHMM::get_posteriors(double**** post)
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	for(int c1Nmod=0; c1Nmod<this->Nmod ; c1Nmod++)
	{
		for(int c2Nmod=0; c2Nmod<this->Nmod ; c2Nmod++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				for (int t=0; t<this->T; t++)
				{
					post[c1Nmod][c2Nmod][iN][t] = this->gamma[c1Nmod][iN][t];
				}
			}
		}
	}

}

//LUISA-changed
double InfluenceScaleHMM::get_posterior(int c1, int iN, int t)
{
	//FILE_LOG(logDEBUG4) << __PRETTY_FUNCTION__;
	return(this->gamma[c1][iN][t]);
}

//--LUISA------solved
double InfluenceScaleHMM::get_density(int c1,int iN, int t)
{
	//FILE_LOG(logDEBUG3) << __PRETTY_FUNCTION__;
	return(this->densities[c1][iN][t]);
}

double InfluenceScaleHMM::get_proba(int c1, int i)
{
	return( this->proba[c1][i] );
}
//--LUISA fixed
double InfluenceScaleHMM::get_A(int c1, int c2, int i, int j)
{
	return( this->A[c1][c2][i][j] );
}
//-----------------
double InfluenceScaleHMM::get_logP()
{
	return( this->logP );
}

void InfluenceScaleHMM::set_cutoff(int cutoff)
{
	this->cutoff = cutoff;
}

// Private ====================================================
// Methods ----------------------------------------------------
void InfluenceScaleHMM::forward()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

//LUISA HERE need to change alpha dimensions !!

		//std::vector<double> alpha(this->N);
		//LUISA -----------------Changed alpha dimension, is it ok? or matrix- not tested yet
		std::vector< std::vector< double > > alpha( this->Nmod, std::vector< double >( this->N ) );
		// Initialization
		//--LUISA-- question :  is alpha initialization right? same chain for densities?

		for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++){

			this->scalefactoralpha[c1Nmod][0] = 0.0;
			for (int iN=0; iN<this->N; iN++)
			{
				alpha[c1Nmod][iN] = this->proba[c1Nmod][iN] * this->densities[c1Nmod][iN][0];
				//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
				this->scalefactoralpha[c1Nmod][0] += alpha[c1Nmod][iN];
				//Rprintf("alpha[c1Nmod=%d][iN=%d](%g) = this->proba[c1Nmod][iN](%g) * this->densities[c1Nmod][iN][0](%g);\n", c1Nmod, iN,alpha[c1Nmod][iN] ,this->proba[c1Nmod][iN] , this->densities[c1Nmod][iN][0] );
			}
		}


		//--LUISA ???
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<0<<"] = " << scalefactoralpha[0];
		for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalealpha[c1Nmod][0][iN] = alpha[c1Nmod][iN] / this->scalefactoralpha[c1Nmod][0];
				//FILE_LOG(logDEBUG4) << "scalealpha["<<0<<"]["<<iN<<"] = " << scalealpha[0][iN];
				//Rprintf("this->scalealpha[c1Nmod=%d][0][iN=%d](%g) = alpha[c1Nmod=%d][iN=%d](%g) / this->scalefactoralpha[c1Nmod=%d][0](%g)\n", c1Nmod, iN,scalealpha[c1Nmod][0][iN],c1Nmod, iN ,alpha[c1Nmod][iN], c1Nmod, scalefactoralpha[c1Nmod][0] );
			}
	  }
		// Induction

				for (int t=1; t<this->T; t++)
				{
					for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
					{
						this->scalefactoralpha[c1Nmod][t] = 0.0;
						for (int iN=0; iN<this->N; iN++)
						{
							double helpsum = 0.0;
							for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
							{
								for (int jN=0; jN<this->N; jN++)
								{
									helpsum += this->scalealpha[c2Nmod][t-1][jN] * this->influence[c2Nmod][c1Nmod][jN][iN];
								}
							}
							alpha[c1Nmod][iN] = helpsum * this->densities[c1Nmod][iN][t];
							//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
							this->scalefactoralpha[c1Nmod][t] += alpha[c1Nmod][iN];

							//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
						}
					}
					//----------
					for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
					{
						for (int iN=0; iN<this->N; iN++)
						{
								this->scalealpha[c1Nmod][t][iN] = alpha[c1Nmod][iN] / this->scalefactoralpha[c1Nmod][t];
								//FILE_LOG(logDEBUG4) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
								if(std::isnan(this->scalealpha[c1Nmod][t][iN]))
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

void InfluenceScaleHMM::backward()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

//LUISA here need to change beta dimensions!!
	std::vector< std::vector< double > > beta( this->Nmod, std::vector< double >( this->N ) );

		// Initialization -LUISA-ok
		for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				beta[c1Nmod][iN] = 1.0;
				//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
			}
	  }
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<T-1<<"] = " << scalefactoralpha[T-1];
		for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalebeta[c1Nmod][T-1][iN] = beta[c1Nmod][iN] / this->scalefactoralpha[c1Nmod][T-1];
				//FILE_LOG(logDEBUG4) << "scalebeta["<<T-1<<"]["<<iN<<"] = " << scalebeta[T-1][iN];
			}
	  }
		// Induction
		for (int t=this->T-2; t>=0; t--)
		{


			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
			{
				for (int iN=0; iN<this->N; iN++)
				{
					beta[c1Nmod][iN] = 0.0;
					for(int c2Nmod=0; c2Nmod<this->Nmod; c2Nmod++)
					{
						for(int jN=0; jN<this->N; jN++)
						{
							//LUISA embedding new formula
							beta[c1Nmod][iN] += this->influence[c1Nmod][c2Nmod][iN][jN]* this->densities[c2Nmod][jN][t+1]*scalebeta[c2Nmod][t+1][jN];
						}
						//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
					}
				}
			}



					//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
			for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
				{
					for (int iN=0; iN<this->N; iN++)
					{
						this->scalebeta[c1Nmod][t][iN] = beta[c1Nmod][iN] / this->scalefactoralpha[c1Nmod][t];
						//FILE_LOG(logDEBUG4) << "scalebeta["<<t<<"]["<<iN<<"] = " << scalebeta[t][iN];
						if (std::isnan(this->scalebeta[c1Nmod][t][iN]))
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

void InfluenceScaleHMM::calc_sumgamma()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	clock_t time = clock(), dtime;

	// // Initialize the sumgamma
	// for(int c1Nmod=0; c1Nmod<this->Nmod; Nmod++)
	// {
	// 		for (int iN=0; iN<this->N; iN++)
	// 		{
	// 			this->sumgamma[c1Nmod][iN] = 0.0;
	// 		}
	// }




	// Compute the gammas (posteriors) and sumgamma
	#pragma omp parallel for
	for(int c1Nmod=0; c1Nmod<this->Nmod; c1Nmod++)
	{
			for (int iN=0; iN<this->N; iN++)
			{
				for (int t=0; t<this->T; t++)
				{
					//--LUISA scalefactor alpha ????

					this->gamma[c1Nmod][iN][t] = this->scalealpha[c1Nmod][t][iN] * this->scalebeta[c1Nmod][t][iN] * this->scalefactoralpha[c1Nmod][t];
				//	this->sumgamma[c1Nmod][iN] += this->gamma[c1Nmod][iN][t];

				}
			}
	}

}
			// Subtract the last value because sumgamma goes only until T-1 and we computed until T to get also loggamma at T
		// 	//--LUISA now for each chain
		// 	for(int c1Nmod=0; c1Nmod<this->Nmod; Nmod++)
		// 	{
		// 		for(int c2Nmod=0; c2Nmod<this->Nmod;Nmod++)
		// 		{
		// 			for (int iN=0; iN<this->N; iN++)
		// 			{
		//
		// 				this->sumgamma[c1Nmod][c2Nmod][iN] -= this->gamma[c1Nmod][c2Nmod][iN][T-1];
		// 			}
		// 		}
		// 	}
		// }


// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_sumgamma(): " << dtime << " clicks";


//------------------LUISA-------------------------------------
void InfluenceScaleHMM::calc_sumxi()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	clock_t time = clock(), dtime;

	double xi;
	// Initialize the sumxi
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
	  {
	  for(int c2Nmod=0; c2Nmod <this->Nmod; c2Nmod++)
	    {
	    for (int iN=0; iN<this->N; iN++)
	    {
	      for (int jN=0; jN<this->N; jN++)
	      {
	        this->sumxi[c1Nmod][c2Nmod][iN][jN] = 0.0;
	      }
	    }

	  }
	}


// 	if (this->xvariate==UNIVARIATE)
// 	{

		#pragma omp parallel for


for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
{
  for(int c2Nmod=0; c2Nmod <this->Nmod; c2Nmod++)
  {
    for (int iN=0; iN<this->N; iN++)
    {
      //FILE_LOG(logDEBUG3) << "Calculating sumxi["<<iN<<"][jN]";
      for (int t=0; t<this->T-1; t++)
      {
        for (int jN=0; jN<this->N; jN++)
        {

          xi = this->scalealpha[c1Nmod][t][iN] * this->A[c1Nmod][c2Nmod][iN][jN] * this->densities[c2Nmod][jN][t+1] * this->scalebeta[c2Nmod][t+1][jN];
          this->sumxi[c1Nmod][c2Nmod][iN][jN] += xi;

        }
      }
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


void InfluenceScaleHMM::calc_tiesumxi(){

	// Initialize
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
		{
		for(int c2Nmod=0; c2Nmod <this->Nmod; c2Nmod++)
			{
				for (int iN=0; iN<this->N; iN++)
				{
					this->tiesumxi[c1Nmod][c2Nmod] = 0.0;
					this->tiesumxitotal[c1Nmod] = 0.0;
					this->tie_onestatus_sumxi[c1Nmod][c2Nmod][iN]=0.0;
				}
		}
	}

	//calculate tiesumxi
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
	{
	  for(int c2Nmod=0; c2Nmod <this->Nmod; c2Nmod++)
	  {
			for (int iN=0; iN<this->N; iN++)
			{
					for (int jN=0; jN<this->N; jN++)
					{
						this->tiesumxi[c1Nmod][c2Nmod]+=this->sumxi[c1Nmod][c2Nmod][iN][jN];
						this->tie_onestatus_sumxi[c1Nmod][c2Nmod][iN]+=this->sumxi[c1Nmod][c2Nmod][iN][jN];
					}
			}
		}
	}

	// calculate tiesumxitotal
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
	{
	  for(int c2Nmod=0; c2Nmod <this->Nmod; c2Nmod++)
	  {
			this->tiesumxitotal[c1Nmod]+= this->tiesumxi[c1Nmod][c2Nmod];
		}
	}

}

//--LUISA here maybe
void InfluenceScaleHMM::calc_loglikelihood()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}// 	clock_t time = clock(), dtime;

	this->logP = 0.0;
	for (int t=0; t<this->T; t++)
	{
		for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
		{
			//only for dimensions!!!!
			this->logP += log(this->scalefactoralpha[c1Nmod][t]);
		}
	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_loglikelihood(): " << dtime << " clicks";
}

void InfluenceScaleHMM::calc_densities()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}//	clock_t time = clock(), dtime;
	// Errors thrown inside a #pragma must be handled inside the thread
	std::vector<bool> nan_encountered(this->N);
	#pragma omp parallel for
	//LUISA : here i changed + only for univariate
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			//FILE_LOG(logDEBUG3) << "Calculating densities for state " << iN;
			try
			{
				//HEREEEE only for dimension, need to check what fr dimension
				//TODO
				//LUISA LUISA LUISA HEY here commented bt not dure if legitim
				this->densityFunctions[c1Nmod][iN]->calc_densities(this->densities[c1Nmod][iN]);
			}
			catch(std::exception& e)
			{
				if (strcmp(e.what(),"nan detected")==0) { nan_encountered[iN]=true; }
				else { throw; }
			}
		}
	}
	for (int iN=0; iN<this->N; iN++)
	{
		if (nan_encountered[iN]==true)
		{
			throw nan_detected;
		}
	}

// ----------------------------------LUISA ------ from here nix mehr
	// Check if the density for all states is close to zero and correct to prevent NaNs
// 	double zero_cutoff = 1.18e-37; // 32-bit precision is 1.18e-38
	double zero_cutoff = 2.23e-307; // 64-bit precision is 2.23e-308
	std::vector<double> temp(this->N);
	// t=0
	for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
	{
	for (int iN=0; iN<this->N; iN++)
	{
		temp[iN] = this->densities[c1Nmod][iN][0];
	}
	}
	if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
	{
		for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
		{
		for (int iN=0; iN<this->N; iN++)
		{
			//here dimensions !!
			this->densities[c1Nmod][iN][0] = zero_cutoff;
		}
		}
	}
	// t>0
	for (int t=1; t<this->T; t++)
	{
		for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
		{
		for (int iN=0; iN<this->N; iN++)
		{
			//HERE
			temp[iN] = this->densities[c1Nmod][iN][t];
		}
		}
		if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
		{
			for(int c1Nmod=0; c1Nmod< this->Nmod; c1Nmod++)
			{
			for (int iN=0; iN<this->N; iN++)
			{
				this->densities[c1Nmod][iN][t] = this->densities[c1Nmod][iN][t-1];
			}
			}
		}
	}

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG) << "calc_densities(): " << dtime << " clicks";
}

/*
void InfluenceScaleHMM::print_uni_iteration(int iteration)
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	if (this->verbosity>=1)
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

*/
void InfluenceScaleHMM::print_multi_iteration(int iteration)
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	if (this->verbosity>=1)
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

/*
void InfluenceScaleHMM::print_uni_params()
{
	if(this->verbosity>=2){ Rprintf("%s\n", __PRETTY_FUNCTION__);}	if (this->verbosity>=2)
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
	*/
//}
