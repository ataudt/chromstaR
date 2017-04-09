#include "R_interface.h"

static ScaleHMM* hmm; // declare as static outside the function because we only need one and this enables memory-cleanup on R_CheckUserInterrupt()
static int** multiO;

// ===================================================================================================================================================
// This function takes parameters from R, creates a univariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// ===================================================================================================================================================
void univariate_hmm(int* O, int* T, int* N, double* size, double* prob, int* maxiter, int* maxtime, double* eps, double* posteriors, double* densities, bool* keep_densities, double* A, double* proba, double* loglik, double* weights, int* iniproc, double* initial_size, double* initial_prob, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff, int* verbosity)
{

	// Define logging level
	//FILE* pFile = fopen("chromStar.log", "w");
 	//Output2FILE::Stream() = pFile;
 	//FILELog::ReportingLevel() = FILELog::FromString("ERROR");
 	//FILELog::ReportingLevel() = FILELog::FromString("DEBUG2");

	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	// Parallelization settings
	#ifdef _OPENMP
	omp_set_num_threads(*num_threads);
	#endif

	// Print some information
	//FILE_LOG(logINFO) << "number of states = " << *N;
	if (*verbosity>=1) Rprintf("HMM: number of states = %d\n", *N);
	//FILE_LOG(logINFO) << "number of bins = " << *T;
	if (*verbosity>=1) Rprintf("HMM: number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		//FILE_LOG(logINFO) << "maximum number of iterations = none";
		if (*verbosity>=1) Rprintf("HMM: maximum number of iterations = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		if (*verbosity>=1) Rprintf("HMM: maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		//FILE_LOG(logINFO) << "maximum running time = none";
		if (*verbosity>=1) Rprintf("HMM: maximum running time = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		if (*verbosity>=1) Rprintf("HMM: maximum running time = %d sec\n", *maxtime);
	}
	//FILE_LOG(logINFO) << "epsilon = " << *eps;
	if (*verbosity>=1) Rprintf("HMM: epsilon = %g\n", *eps);

	//FILE_LOG(logDEBUG3) << "observation vector";
	for (int t=0; t<50; t++) {
		//FILE_LOG(logDEBUG3) << "O["<<t<<"] = " << O[t];
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Create the HMM
	//FILE_LOG(logDEBUG1) << "Creating a univariate HMM";
	hmm = new ScaleHMM(*T, *N, *verbosity);
	hmm->set_cutoff(*read_cutoff);
	// Initialize the transition probabilities and proba
	hmm->initialize_transition_probs(initial_A, *use_initial_params);
	hmm->initialize_proba(initial_proba, *use_initial_params);
    
	// Calculate mean and variance of data
	double Tadjust = 0, mean = 0, variance = 0;
	for(int t=0; t<*T; t++)
	{
		if (O[t]>0)
		{
			mean += O[t];
			Tadjust += 1;
		}
	}
	mean = mean / Tadjust;
	for(int t=0; t<*T; t++)
	{
		if (O[t]>0)
		{
			variance += pow(O[t] - mean, 2);
		}
	}
	variance = variance / Tadjust;
	//FILE_LOG(logINFO) << "data mean = " << mean << ", data variance = " << variance;		
	if (*verbosity>=1) Rprintf("HMM: data mean = %g, data variance = %g\n", mean, variance);		
	
	// Go through all states of the hmm and assign the density functions
	double imean=0, ivariance=0;
	for (int i_state=0; i_state<*N; i_state++)
	{
		if (*use_initial_params) {
			//FILE_LOG(logINFO) << "Using given parameters for size and prob";
			if (*verbosity>=1) Rprintf("HMM: Using given parameters for size and prob\n");
			imean = (1-initial_prob[i_state])*initial_size[i_state] / initial_prob[i_state];
			ivariance = imean / initial_prob[i_state];
			//FILE_LOG(logDEBUG2) << "imean = " << imean;
			//FILE_LOG(logDEBUG2) << "ivariance = " << ivariance;
		} else {

			if (*iniproc == 1)
			{
				// Simple initialization, seems to give the fastest convergence
				if (i_state == 1)
				{
					//FILE_LOG(logDEBUG) << "Initializing size and prob for state 1";
					imean = mean;
					ivariance = variance;
				}
				else if (i_state == 2)
				{
					//FILE_LOG(logDEBUG) << "Initializing size and prob for state 2";
					imean = mean*1.5;
					ivariance = variance*2;
				} 
				// Make sure variance is greater than mean
				if (imean >= ivariance)
				{
					ivariance = imean + 1;
				}
			}
			else if (*iniproc == 2)
			{
				// Disturb mean and variance for use as randomized initial parameters
				//FILE_LOG(logINFO) << "Using random initialization for size and prob";
				if (*verbosity>=1) Rprintf("HMM: Using random initialization for size and prob\n");
					imean = runif(0, 10*mean);
					ivariance = imean + runif(0, 20*imean); // variance has to be greater than mean, otherwise r will be negative
				//FILE_LOG(logDEBUG2) << "imean = " << imean;
				//FILE_LOG(logDEBUG2) << "ivariance = " << ivariance;
			}
			else if (*iniproc == 3)
			{
				// Empirical initialization
				if (i_state == 1)
				{
					//FILE_LOG(logINFO) << "Initializing r and p empirically for state 1";
					if (*verbosity>=1) Rprintf("HMM: Initializing r and p empirically for state 1\n");
					imean = mean/2;
					ivariance = imean*2;
				}
				else if (i_state == 2)
				{
					//FILE_LOG(logINFO) << "Initializing r and p empirically for state 2";
					if (*verbosity>=1) Rprintf("HMM: Initializing r and p empirically for state 2\n");
					imean = mean*2;
					ivariance = imean*2;
				} 
			}

			// Calculate r and p from mean and variance
			initial_size[i_state] = pow(imean,2)/(ivariance-imean);
			initial_prob[i_state] = imean/ivariance;

		}

		if (i_state >= 1)
		{
			//FILE_LOG(logDEBUG1) << "Using negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, *T, initial_size[i_state], initial_prob[i_state]); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else if (i_state == 0)
		{
			//FILE_LOG(logDEBUG1) << "Using only zeros for state " << i_state;
			ZeroInflation *d = new ZeroInflation(O, *T); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else
		{
			//FILE_LOG(logWARNING) << "Density not specified, using default negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, *T, initial_size[i_state], initial_prob[i_state]);
			hmm->densityFunctions.push_back(d);
		}
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Do the Baum-Welch to estimate the parameters
	//FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		hmm->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		//FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		if (*verbosity>=1) Rprintf("HMM: Error in Baum-Welch: %s\n", e.what());
		if (strcmp(e.what(),"nan detected")==0) { *error = 1; }
		else { *error = 2; }
	}
		
	//FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";

	// Get the posteriors and save results directly to the R pointer
	//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
	if (*verbosity>=1) Rprintf("HMM: Recoding posteriors ...\n");
	R_FlushConsole();
	#pragma omp parallel for
	for (int iN=0; iN<*N; iN++)
	{
		for (int t=0; t<*T; t++)
		{
			posteriors[t + iN * (*T)] = hmm->get_posterior(iN, t);
		}
	}

	// Get the densities and save results directly to the R pointer
	if (*keep_densities == true)
	{
		//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
		if (*verbosity>=1) Rprintf("HMM: Recoding densities ...\n");
		R_FlushConsole();
		#pragma omp parallel for
		for (int iN=0; iN<*N; iN++)
		{
			for (int t=0; t<*T; t++)
			{
				densities[t + iN * (*T)] = hmm->get_density(iN, t);
			}
		}
	}

	//FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<*N; i++)
	{
		proba[i] = hmm->get_proba(i);
		for (int j=0; j<*N; j++)
		{
			A[i * (*N) + j] = hmm->get_A(j,i);
		}
	}

	// copy the estimated distribution params
	for (int i=0; i<*N; i++)
	{
		if (hmm->densityFunctions[i]->get_name() == NEGATIVE_BINOMIAL) 
		{
			NegativeBinomial* d = (NegativeBinomial*)(hmm->densityFunctions[i]);
			size[i] = d->get_size();
			prob[i] = d->get_prob();
		}
		else if (hmm->densityFunctions[i]->get_name() == ZERO_INFLATION)
		{
			// These values for a Negative Binomial define a zero-inflation (delta distribution)
			size[i] = 0;
			prob[i] = 1;
		}
	}
	*loglik = hmm->get_logP();
	hmm->calc_weights(weights);
	
	//FILE_LOG(logDEBUG1) << "Deleting the hmm";
	delete hmm;
	hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call
}

// =====================================================================================================================================================
// This function takes parameters from R, creates a multivariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// =====================================================================================================================================================
void multivariate_hmm(int* O, int* T, int* N, int *Nmod, double* comb_states, double* size, double* prob, double* w, double* cor_matrix_inv, double* det, int* maxiter, int* maxtime, double* eps, double* posteriors, bool* keep_posteriors, double* densities, bool* keep_densities, int* states, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* verbosity)
{

	// Define logging level {"ERROR", "WARNING", "INFO", "ITERATION", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"}
 	//FILE* pFile = fopen("chromStar.log", "w");
	//Output2FILE::Stream() = pFile;
 	//FILELog::ReportingLevel() = FILELog::FromString("ERROR");
 	//FILELog::ReportingLevel() = FILELog::FromString("DEBUG3");

	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	// Parallelization settings
	#ifdef _OPENMP
	omp_set_num_threads(*num_threads);
	#endif

	// Print some information
	//FILE_LOG(logINFO) << "number of states = " << *N;
	if (*verbosity>=1) Rprintf("HMM: number of states = %d\n", *N);
	//FILE_LOG(logINFO) << "number of bins = " << *T;
	if (*verbosity>=1) Rprintf("HMM: number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		//FILE_LOG(logINFO) << "maximum number of iterations = none";
		if (*verbosity>=1) Rprintf("HMM: maximum number of iterations = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		if (*verbosity>=1) Rprintf("HMM: maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		//FILE_LOG(logINFO) << "maximum running time = none";
		if (*verbosity>=1) Rprintf("HMM: maximum running time = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		if (*verbosity>=1) Rprintf("HMM: maximum running time = %d sec\n", *maxtime);
	}
	//FILE_LOG(logINFO) << "epsilon = " << *eps;
	if (*verbosity>=1) Rprintf("HMM: epsilon = %g\n", *eps);
	//FILE_LOG(logINFO) << "number of experiments = " << *Nmod;
	if (*verbosity>=1) Rprintf("HMM: number of experiments = %d\n", *Nmod);

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Recode the observation vector to matrix representation
// 	clock_t clocktime = clock(), dtime;
	multiO = CallocIntMatrix(*Nmod, *T);
	for (int imod=0; imod<*Nmod; imod++)
	{
		for (int t=0; t<*T; t++)
		{
			multiO[imod][t] = O[imod*(*T)+t];
		}
	}
// 	dtime = clock() - clocktime;
// 	//FILE_LOG(logDEBUG1) << "recoding observation vector to matrix representation: " << dtime << " clicks";

	// Create the HMM
	//FILE_LOG(logDEBUG1) << "Creating the multivariate HMM";
	hmm = new ScaleHMM(*T, *N, *Nmod, *verbosity);
	// Initialize the transition probabilities and proba
	hmm->initialize_transition_probs(initial_A, *use_initial_params);
	hmm->initialize_proba(initial_proba, *use_initial_params);
	
	// Print logproba and A
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		//FILE_LOG(logDEBUG) << "proba["<<iN<<"] = " <<exp(hmm->logproba[iN]);
// 		for (int jN=0; jN<*N; jN++)
// 		{
// 			//FILE_LOG(logDEBUG) << "A["<<iN<<"]["<<jN<<"] = " << hmm->A[iN][jN];
// 		}
// 	}

	// Prepare the binary_states (univariate) vector: binary_states[N][Nmod], e.g., binary_states[iN][imod] tells me at state comb_states[iN], modification imod is non-enriched (0) or enriched (1)
	//FILE_LOG(logDEBUG1) << "Preparing the binary_states vector";
	double res;
	bool **binary_states = CallocBoolMatrix(*N, *Nmod);
	for (int iN=0; iN < *N; iN++) //for each comb state considered
	{
		res = comb_states[iN];
		for (int imod=(*Nmod-1); imod >= 0; imod--) //for each modification of this comb state
		{
			binary_states[iN][imod] = (bool)fmod(res,2);
			res = (res - (double)binary_states[iN][imod]) / 2.0;
		}
	}

	/* initialize the distributions */
	//FILE_LOG(logDEBUG1) << "Initializing the distributions";
	for (int iN=0; iN<*N; iN++) //for each combinatorial state
	{
		std::vector <Density*> tempMarginals;            
		for (int imod=0; imod < *Nmod; imod++) //for each modification
		{
			Density *d;
			if (binary_states[iN][imod]) //construct the marginal density function for modification imod being enriched
			{
				d = new NegativeBinomial(multiO[imod], *T, size[2*imod+1], prob[2*imod+1]); // delete is done inside ~MVCopulaApproximation()
			}
			else //construct the density function for modification imod being non-enriched
			{
				d = new ZiNB(multiO[imod], *T, size[2*imod], prob[2*imod], w[imod]); // delete is done inside ~MVCopulaApproximation()
			}
			tempMarginals.push_back(d);
		}
		//MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(O, tempMarginals, &(cor_matrix_inv[iN*Nmod*Nmod]), det[iN]);
		//FILE_LOG(logDEBUG1) << "Calling MVCopulaApproximation for state " << iN;
		MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(multiO, *T, tempMarginals, &(cor_matrix_inv[iN**Nmod**Nmod]), det[iN]); // delete is done inside ~ScaleHMM()
		hmm->densityFunctions.push_back(tempMVdens);
	}
	FreeBoolMatrix(binary_states, *N);
	
	// Estimate the parameters
	//FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		hmm->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		//FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		if (*verbosity>=1) Rprintf("HMM: Error in Baum-Welch: %s\n", e.what());
		if (strcmp(e.what(),"nan detected")==0) { *error = 1; }
		else { *error = 2; }
	}
	//FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	
	// Get the posteriors and save results directly to the R pointer
	if (*keep_posteriors == true)
	{
		//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
		if (*verbosity>=1) Rprintf("HMM: Recoding posteriors ...\n");
		R_FlushConsole();
		#pragma omp parallel for
		for (int iN=0; iN<*N; iN++)
		{
			for (int t=0; t<*T; t++)
			{
				posteriors[t + iN * (*T)] = hmm->get_posterior(iN, t);
			}
		}
	}

	// Get the densities and save results directly to the R pointer
	if (*keep_densities == true)
	{
		//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
		if (*verbosity>=1) Rprintf("HMM: Recoding densities ...\n");
		R_FlushConsole();
		#pragma omp parallel for
		for (int iN=0; iN<*N; iN++)
		{
			for (int t=0; t<*T; t++)
			{
				densities[t + iN * (*T)] = hmm->get_density(iN, t);
			}
		}
	}

	// Compute the states from posteriors
	//FILE_LOG(logDEBUG1) << "Computing states from posteriors";
// 	if (*fdr == -1)
// 	{
		int ind_max;
		std::vector<double> posterior_per_t(*N);
		for (int t=0; t<*T; t++)
		{
			for (int iN=0; iN<*N; iN++)
			{
				posterior_per_t[iN] = hmm->get_posterior(iN, t);
			}
			ind_max = std::distance(posterior_per_t.begin(), std::max_element(posterior_per_t.begin(), posterior_per_t.end()));
			states[t] = comb_states[ind_max];
		}
// 	}
// 	else
// 	{
// 		double** transformed_posteriors = CallocDoubleMatrix(*T, *Nmod);
// 		for (int t=0; t<*T; t++)
// 		{
// 			for (int iN=0; iN<*N; iN++)
// 			{
// 				for (int iNmod=0; iNmod<*Nmod; iNmod++)
// 				{
// 					transformed_posteriors[t][iNmod] += (double)binary_states[iN][iNmod] * hmm->get_posterior(iN, t);
// 				}
// 			}
// 		}
// 	}
	
	//FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<*N; i++)
	{
		proba[i] = hmm->get_proba(i);
		for (int j=0; j<*N; j++)
		{
				A[i * (*N) + j] = hmm->get_A(i,j);
		}
	}
	*loglik = hmm->get_logP();

	//FILE_LOG(logDEBUG1) << "Deleting the hmm";
	delete hmm;
	hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call
// 	FreeIntMatrix(multiO, *Nmod); // free on.exit() in R code
}


// =======================================================
// This function make a cleanup if anything was left over
// =======================================================
void univariate_cleanup()
{
// 	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__; // This message will be shown if interrupt happens before start of C-code
	delete hmm;
}

void multivariate_cleanup(int* Nmod)
{
	delete hmm;
	FreeIntMatrix(multiO, *Nmod);
}


// =======================================================
// C version of apply(array3D, 1, which.max)
// =======================================================
void array3D_which_max(double* array3D, int* dim, int* ind_max)
{
  // array3D is actually a vector, but is intended to originate from a 3D array in R
	std::vector<double> value_per_i0(dim[1] * dim[2]);
  for (int i0=0; i0<dim[0]; i0++)
  {
    for (int i1=0; i1<dim[1]; i1++)
    {
      for (int i2=0; i2<dim[2]; i2++)
      {
  			value_per_i0[i1*dim[2]+i2] = array3D[(i1*dim[2]+i2) * dim[0] + i0];
        // Rprintf("i0=%d, i1=%d, i2=%d, value_per_i0[%d] = %g\n", i0, i1, i2, i1*dim[2]+i2, value_per_i0[i1*dim[2]+i2]);
      }
    }
		ind_max[i0] = 1 + std::distance(value_per_i0.begin(), std::max_element(value_per_i0.begin(), value_per_i0.end()));
  }
	
}


// ====================================================================================
// C version of apply(array2D, MARGIN = 1, FUN = mean)
// ====================================================================================
void array2D_mean(double* array2D, int* dim, double* mean)
{
  // array2D is actually a vector, but is intended to originate from a matrix in R
  double sum=0;
  for (int i0=0; i0<dim[0]; i0++)
  {
    sum = 0;
    for (int i1=0; i1<dim[1]; i1++)
    {
      sum += array2D[i1*dim[0] + i0];
    }
    mean[i0] = sum / dim[1];
  }
	
}


// ====================================================================================
// C version of apply(array3D, MARGIN = c(1,2), FUN = mean)
// ====================================================================================
void array3D_mean(double* array3D, int* dim, double* mean)
{
  // array3D is actually a vector, but is intended to originate from a 3D array in R
  double sum=0;
  for (int i0=0; i0<dim[0]; i0++)
  {
    for (int i1=0; i1<dim[1]; i1++)
    {
      sum = 0;
      for (int i2=0; i2<dim[2]; i2++)
      {
        sum += array3D[(i2*dim[1] + i1)*dim[0] + i0];
        // Rprintf("i0=%d, i1=%d, i2=%d, array3D[(i2*dim[1] + i1)*dim[0] + i0 = %d] = %g\n", i0, i1, i2, (i2*dim[1] + i1)*dim[0] + i0, array3D[(i2*dim[1] + i1)*dim[0] + i0]);
      }
      mean[i1*dim[0] + i0] = sum / dim[2];
    }
  }
	
}