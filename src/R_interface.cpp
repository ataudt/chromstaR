#include "utility.h"
#include "scalehmm.h"
#include <vector> // storing density functions in multivariate
#include <omp.h> // parallelization options

// ===================================================================================================================================================
// This function takes parameters from R, creates a univariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// ===================================================================================================================================================
extern "C" {
void R_univariate_hmm(int* O, int* T, int* N, double* r, double* p, int* maxiter, int* maxtime, double* eps, double* posteriors, double* A, double* proba, double* loglik, double* weights, int* iniproc, double* initial_r, double* initial_p, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff)
{

	// Define logging level
// 	FILE* pFile = fopen("chromStar.log", "w");
// 	Output2FILE::Stream() = pFile;
 	FILELog::ReportingLevel() = FILELog::FromString("ERROR");

	// Parallelization settings
	omp_set_num_threads(*num_threads);

	// Print some information
	FILE_LOG(logINFO) << "number of states = " << *N;
	Rprintf("number of states = %d\n", *N);
	FILE_LOG(logINFO) << "number of bins = " << *T;
	Rprintf("number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		FILE_LOG(logINFO) << "maximum number of iterations = none";
		Rprintf("maximum number of iterations = none\n");
	} else {
		FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		Rprintf("maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		FILE_LOG(logINFO) << "maximum running time = none";
		Rprintf("maximum running time = none\n");
	} else {
		FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		Rprintf("maximum running time = %d sec\n", *maxtime);
	}
	FILE_LOG(logINFO) << "epsilon = " << *eps;
	Rprintf("epsilon = %g\n", *eps);

	FILE_LOG(logDEBUG3) << "observation vector";
	for (int t=0; t<50; t++) {
		FILE_LOG(logDEBUG3) << "O["<<t<<"] = " << O[t];
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Create the HMM
	FILE_LOG(logDEBUG1) << "Creating a univariate HMM";
	ScaleHMM* model = new ScaleHMM(*T, *N);
	model->set_cutoff(*read_cutoff);
	// Initialize the transition probabilities and proba
	model->initialize_transition_probs(initial_A, *use_initial_params);
	model->initialize_proba(initial_proba, *use_initial_params);
    
	// Calculate mean and variance of data
	double mean = 0, variance = 0;
	for(int t=0; t<*T; t++)
	{
		mean+= O[t];
	}
	mean = mean / *T;
	for(int t=0; t<*T; t++)
	{
		variance+= pow(O[t] - mean, 2);
	}
	variance = variance / *T;
	FILE_LOG(logINFO) << "data mean = " << mean << ", data variance = " << variance;		
	Rprintf("data mean = %g, data variance = %g\n", mean, variance);		
	
	// Go through all states of the model and assign the density functions
	double imean, ivariance;
	for (int i_state=0; i_state<model->N; i_state++)
	{
		if (*use_initial_params) {
			FILE_LOG(logINFO) << "Using given parameters for r and p";
			Rprintf("Using given parameters for r and p\n");
			imean = (1-initial_p[i_state])*initial_r[i_state] / initial_p[i_state];
			ivariance = imean / initial_p[i_state];
			FILE_LOG(logDEBUG2) << "imean = " << imean;
			FILE_LOG(logDEBUG2) << "ivariance = " << ivariance;
		} else {

			if (*iniproc == 1)
			{
				// Simple initialization, seems to give the fastest convergence
				if (i_state == 1)
				{
					FILE_LOG(logDEBUG) << "Initializing r and p for state 1";
					imean = mean;
					ivariance = variance;
				}
				else if (i_state == 2)
				{
					FILE_LOG(logDEBUG) << "Initializing r and p for state 2";
					imean = mean+1;
					ivariance = variance;
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
				FILE_LOG(logINFO) << "Using random initialization for r and p";
				Rprintf("Using random initialization for r and p\n");
				srand (clock());
				int rand1, rand2;
				rand1 = rand();
				rand2 = rand();
				imean = (double)rand1/(double)RAND_MAX * 10*mean;
				ivariance = imean + (double)rand2/(double)RAND_MAX * 20*imean; // variance has to be greater than mean, otherwise r will be negative
				FILE_LOG(logDEBUG2) << "RAND_MAX = " << RAND_MAX;
				FILE_LOG(logDEBUG2) << "rand1 = " << rand1;
				FILE_LOG(logDEBUG2) << "rand2 = " << rand2;
				FILE_LOG(logDEBUG2) << "imean = " << imean;
				FILE_LOG(logDEBUG2) << "ivariance = " << ivariance;
			}
			else if (*iniproc == 3)
			{
				// Empirical initialization
				if (i_state == 1)
				{
					FILE_LOG(logINFO) << "Initializing r and p empirically for state 1";
					Rprintf("Initializing r and p empirically for state 1\n");
					imean = mean/2;
					ivariance = imean*2;
				}
				else if (i_state == 2)
				{
					FILE_LOG(logINFO) << "Initializing r and p empirically for state 2";
					Rprintf("Initializing r and p empirically for state 2\n");
					imean = mean*2;
					ivariance = imean*2;
				} 
			}

			// Calculate r and p from mean and variance
			initial_r[i_state] = pow(imean,2)/(ivariance-imean);
			initial_p[i_state] = imean/ivariance;

		}

		if (i_state >= 1)
		{
			FILE_LOG(logDEBUG1) << "Using negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, model->T, initial_r[i_state], initial_p[i_state]); // delete is done inside ~LogHMM()
			model->densityFunctions.push_back(d);
		}
		else if (i_state == 0)
		{
			FILE_LOG(logDEBUG1) << "Using only zeros for state " << i_state;
			OnlyZeros *d = new OnlyZeros(O, model->T); // delete is done inside ~LogHMM()
			model->densityFunctions.push_back(d);
		}
		else
		{
			FILE_LOG(logWARNING) << "Density not specified, using default negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, model->T, initial_r[i_state], initial_p[i_state]);
			model->densityFunctions.push_back(d);
		}
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Do the Baum-Welch to estimate the parameters
	FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		model->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		Rprintf("Error in Baum-Welch: %s\n", e.what());
		if (e.what()=="nan detected") { *error = 1; }
		else { *error = 2; }
	}
		
	FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	// Compute the posteriors and save results directly to the R pointer
	double** post = allocDoubleMatrix(model->N, model->T);
	model->get_posteriors(post);
	FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
	#pragma omp parallel for
	for (int iN=0; iN<model->N; iN++)
	{
		for (int t=0; t<model->T; t++)
		{
			posteriors[t + iN * model->T] = post[iN][t];
		}
	}
	freeDoubleMatrix(post, model->N);

	FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<model->N; i++)
	{
		proba[i] = model->get_proba(i);
		for (int j=0; j<model->N; j++)
		{
          A[i*model->N + j] = model->A[i][j];
		}
	}

	// copy the estimated distribution params
	for (int i=0; i<model->N; i++)
	{
		if (model->densityFunctions[i]->getType() == NB) 
		{
			NegativeBinomial* d = (NegativeBinomial*)(model->densityFunctions[i]);
			r[i] = d->getR();
			p[i] = d->getP();
		}
		else if (model->densityFunctions[i]->getType() == ZI)
		{
			// These values for a Negative Binomial define a zero-inflation (delta distribution)
			r[i] = 0;
			p[i] = 1;
		}
	}
	*loglik = model->logP;
	model->calc_weights(weights);
	
	FILE_LOG(logDEBUG1) << "Deleting the model";
	delete model;
}
} // extern C

// =====================================================================================================================================================
// This function takes parameters from R, creates a multivariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// =====================================================================================================================================================
extern "C" {
void R_multivariate_hmm(int* O, int* T, int* N, int *Nmod, int* states, double* r, double* p, double* w, double* cor_matrix_inv, double* det, int* maxiter, int* maxtime, double* eps, double* posteriors, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error)
{

	// Define logging level {"ERROR", "WARNING", "INFO", "ITERATION", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"}
// 	FILE* pFile = fopen("chromStar.log", "w");
// 	Output2FILE::Stream() = pFile;
//  	FILELog::ReportingLevel() = FILELog::FromString("ITERATION");
//  	FILELog::ReportingLevel() = FILELog::FromString("DEBUG");
 	FILELog::ReportingLevel() = FILELog::FromString("ERROR");

	// Parallelization settings
// 	omp_set_num_threads(*num_threads);

	// Print some information
	FILE_LOG(logINFO) << "number of states = " << *N;
	Rprintf("number of states = %d\n", *N);
	FILE_LOG(logINFO) << "number of bins = " << *T;
	Rprintf("number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		FILE_LOG(logINFO) << "maximum number of iterations = none";
		Rprintf("maximum number of iterations = none\n");
	} else {
		FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		Rprintf("maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		FILE_LOG(logINFO) << "maximum running time = none";
		Rprintf("maximum running time = none\n");
	} else {
		FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		Rprintf("maximum running time = %d sec\n", *maxtime);
	}
	FILE_LOG(logINFO) << "epsilon = " << *eps;
	Rprintf("epsilon = %g\n", *eps);
	FILE_LOG(logINFO) << "number of modifications = " << *Nmod;
	Rprintf("number of modifications = %d\n", *Nmod);

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Recode the observation vector to matrix representation
// 	clock_t clocktime = clock(), dtime;
	int** multiO = allocIntMatrix(*Nmod, *T);
	for (int imod=0; imod<*Nmod; imod++)
	{
		for (int t=0; t<*T; t++)
		{
			multiO[imod][t] = O[imod*(*T)+t];
		}
	}
// 	dtime = clock() - clocktime;
// 	FILE_LOG(logDEBUG1) << "recoding observation vector to matrix representation: " << dtime << " clicks";

	// Create the HMM
	FILE_LOG(logDEBUG1) << "Creating the multivariate HMM";
// 	LogHMM* model = new LogHMM(*T, *N, *Nmod);
	ScaleHMM* model = new ScaleHMM(*T, *N, *Nmod);
	// Initialize the transition probabilities and proba
	model->initialize_transition_probs(initial_A, *use_initial_params);
	model->initialize_proba(initial_proba, *use_initial_params);
	
	// Print logproba and A
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		FILE_LOG(logDEBUG) << "proba["<<iN<<"] = " <<exp(model->logproba[iN]);
// 		for (int jN=0; jN<*N; jN++)
// 		{
// 			FILE_LOG(logDEBUG) << "A["<<iN<<"]["<<jN<<"] = " << model->A[iN][jN];
// 		}
// 	}

	// Prepare the binary_states (univariate) vector: binary_states[N][Nmod], e.g., binary_states[iN][imod] tells me at state states[iN], modification imod is non-enriched (0) or enriched (1)
	FILE_LOG(logDEBUG1) << "Preparing the binary_states vector";
	bool **binary_states = allocBoolMatrix(model->N, model->Nmod);
	for(int iN=0; iN < model->N; iN++) //for each comb state considered
	{
		for(int imod=0; imod < model->Nmod; imod++) //for each modification of this comb state
		{
			binary_states[iN][imod] = states[iN]&(int)pow(2,model->Nmod-imod-1);//if =0->hidden state states[iN] has modification imod non enriched; if !=0->enriched
			if (binary_states[iN][imod] != 0 )
				binary_states[iN][imod] = 1;
		}
	}

	/* initialize the distributions */
	FILE_LOG(logDEBUG1) << "Initializing the distributions";
	for (int iN=0; iN<model->N; iN++) //for each combinatorial state
	{
		std::vector <Density*> tempMarginals;            
		for (int imod=0; imod < model->Nmod; imod++) //for each modification
		{
			Density *d;
			if (binary_states[iN][imod]) //construct the marginal density function for modification imod being enriched
			{
				d = new NegativeBinomial(multiO[imod], model->T, r[2*imod+1], p[2*imod+1]); // delete is done inside ~MVCopulaApproximation()
			}
			else //construct the density function for modification imod being non-enriched
			{
				d = new ZiNB(multiO[imod], model->T, r[2*imod], p[2*imod], w[imod]); // delete is done inside ~MVCopulaApproximation()
			}
			tempMarginals.push_back(d);
		}
		//MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(O, tempMarginals, &(cor_matrix_inv[iN*Nmod*Nmod]), det[iN]);
		FILE_LOG(logDEBUG1) << "Calling MVCopulaApproximation for state " << iN;
		MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(multiO, model->T, tempMarginals, &(cor_matrix_inv[iN*model->Nmod*model->Nmod]), det[iN]); // delete is done inside ~LogHMM()
		model->densityFunctions.push_back(tempMVdens);
	}
	
	// Estimate the parameters
	FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		model->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		Rprintf("Error in Baum-Welch: %s\n", e.what());
		if (e.what()=="nan detected") { *error = 1; }
		else { *error = 2; }
	}
	FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	
	// Compute the posteriors and save results directly to the R pointer
	double** post = allocDoubleMatrix(model->N, model->T);
	model->get_posteriors(post);
	FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
	for (int iN=0; iN<model->N; iN++)
	{
		for (int t=0; t<model->T; t++)
		{
			posteriors[t + iN * model->T] = post[iN][t];
		}
	}
	freeDoubleMatrix(post, model->N);
	
	FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<model->N; i++)
	{
		proba[i] = model->get_proba(i);
		for (int j=0; j<model->N; j++)
		{
				A[i*model->N + j] = model->A[i][j];
		}
	}
	*loglik = model->logP;

	FILE_LOG(logDEBUG1) << "Deleting the model";
	delete model;
	freeIntMatrix(multiO, *Nmod);
	freeBoolMatrix(binary_states, *N);
}
} // extern C

// ---------------------------------------------------------------
// void R_multivariate_hmm_productBernoulli()
// This function takes parameters from R, creates a multivariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// ---------------------------------------------------------------
extern "C" {//observation is now the posterior for being UNMODIFIED
void R_multivariate_hmm_productBernoulli(double* O, int* T, int* N, int *Nmod, int* states, int* maxiter, int* maxtime, double* eps, double* post, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error)
{

	// Define logging level
// 	FILELog::ReportingLevel() = FILELog::FromString("DEBUG3");
	FILELog::ReportingLevel() = FILELog::FromString("ITERATION");

	// Parallelization settings
// 	omp_set_num_threads(*num_threads);

	// Print some information
	FILE_LOG(logINFO) << "number of states = " << *N;
	FILE_LOG(logINFO) << "number of bins = " << *T;
	if (*maxiter < 0)
	{
		FILE_LOG(logINFO) << "maximum number of iterations = none";
	} else {
		FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
	}
	if (*maxtime < 0)
	{
		FILE_LOG(logINFO) << "maximum running time = none";
	} else {
		FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
	}
	FILE_LOG(logINFO) << "epsilon = " << *eps;
	FILE_LOG(logINFO) << "number of modifications = " << *Nmod;

	// Recode the observation vector to matrix representation
// 	clock_t clocktime = clock(), dtime;
	double** multiO = allocDoubleMatrix(*Nmod, *T);
	for (int imod=0; imod<*Nmod; imod++)
	{
		for (int t=0; t<*T; t++)
		{
			multiO[imod][t] = O[imod*(*T)+t];
		}
	}
// 	dtime = clock() - clocktime;
// 	FILE_LOG(logDEBUG1) << "recoding observation and probability vectors to matrix representation: " << dtime << " clicks";

	// Create the HMM
	FILE_LOG(logDEBUG1) << "Creating the multivariate HMM";
// 	LogHMM* model = new LogHMM(*T, *N, *Nmod);
	ScaleHMM* model = new ScaleHMM(*T, *N, *Nmod);

	// Initialize the transition probabilities and proba
	model->initialize_transition_probs(initial_A, *use_initial_params);
	model->initialize_proba(initial_proba, *use_initial_params);
	
	// Print logproba and A
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		FILE_LOG(logINFO) << "proba["<<iN<<"] = " <<exp(model->logproba[iN]);
// 		for (int jN=0; jN<*N; jN++)
// 		{
// 			FILE_LOG(logINFO) << "A["<<iN<<"]["<<jN<<"] = " << model->A[iN][jN];
// 		}
// 	}

	// Prepare the binary_states (univariate) vector: binary_states[N][Nmod], e.g., binary_states[iN][imod] tells me at state states[iN], modification imod is non-enriched (0) or enriched (1)
	FILE_LOG(logDEBUG1) << "Preparing the binary_states vector";
	bool **binary_states = allocBoolMatrix(model->N, model->Nmod);
	for(int iN=0; iN < model->N; iN++) //for each comb state considered
	{
		for(int imod=0; imod < model->Nmod; imod++) //for each modification of this comb state
		{
			binary_states[iN][imod] = states[iN]&(int)pow(2,model->Nmod-imod-1);//if =0->hidden state states[iN] has modification imod non enriched; if !=0->enriched
			if (binary_states[iN][imod] != 0 )
				binary_states[iN][imod] = 1;
		}
	}

	/* initialize the distributions */
	FILE_LOG(logDEBUG1) << "Initializing the distributions";
	for (int iN=0; iN<model->N; iN++) //for each combinatorial state
	{
		FILE_LOG(logDEBUG1) << "Calling BernoulliProduct";
		BernoulliProduct *tempBP = new BernoulliProduct(multiO, binary_states[iN], *T, *Nmod); // delete is done inside ~LogHMM()
		model->densityFunctions.push_back(tempBP);
	}

	// Estimate the parameters
	FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		model->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		Rprintf("Error in Baum-Welch: %s\n", e.what());
		if (e.what()=="nan detected") { *error = 1; }
		else { *error = 2; }
	}
	FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	
	// Compute the posteriors and save results directly to the R pointer
	double** posteriors = allocDoubleMatrix(model->N, model->T);
	model->get_posteriors(posteriors);
	FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
	for (int iN=0; iN<model->N; iN++)
	{
		for (int t=0; t<model->T; t++)
		{
			post[t + iN * model->T] = posteriors[iN][t];
		}
	}
	freeDoubleMatrix(posteriors, model->N);
	
	FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<model->N; i++)
	{
		proba[i] = model->get_proba(i);
		for (int j=0; j<model->N; j++)
		{
				A[i*model->N + j] = model->A[i][j];
		}
	}
	*loglik = model->logP;

	FILE_LOG(logDEBUG1) << "Deleting the model";
	delete model;
	freeDoubleMatrix(multiO, *Nmod);
	freeBoolMatrix(binary_states, *N);
}
} // extern C
