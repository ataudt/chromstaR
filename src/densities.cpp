#include "densities.h"

// ============================================================
// Zero-inflated Negative Binomial
// ============================================================

// Constructor and Destructor ---------------------------------
ZiNB::ZiNB()
{
	this->lxfactorials = NULL;
}

ZiNB::ZiNB(int* observations, int T, double size, double prob, double w)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->prob = prob;
	this->size = size;
	this->w = w;
	this->lxfactorials = NULL;
	if (this->obs != NULL)
	{
		this->max_obs = intMax(observations, T);
		this->lxfactorials = (double*) Calloc(max_obs+1, double);
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_obs; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

ZiNB::~ZiNB()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->lxfactorials != NULL)
	{
		Free(this->lxfactorials);
	}
}

// Methods ----------------------------------------------------
void ZiNB::calc_logdensities(double* logdens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		std::vector<double> lGammaRplusX(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			lGammaRplusX[j] = lgamma(this->size + j);
		}
		for (int t=0; t<this->T; t++)
		{
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				logdens[t] = log( this->w + (1-this->w) * exp( lGammaRplusX[(int) this->obs[t]] - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp ) );
			}
			else
			{
				logdens[t] = log(1-this->w) + lGammaRplusX[(int) this->obs[t]] - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			}
			if (std::isnan(logdens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		//FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				logdens[t] = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp ) );
			}
			else
			{
				logdens[t] = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			}
			if (std::isnan(logdens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_densities(double* dens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		std::vector<double> lGammaRplusX(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			lGammaRplusX[j] = lgamma(this->size + j);
		}
		for (int t=0; t<this->T; t++)
		{
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				dens[t] = this->w + (1-this->w) * exp( lGammaRplusX[(int) this->obs[t]] - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			else
			{
				dens[t] = (1-this->w) * exp( lGammaRplusX[(int) this->obs[t]] - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			if (std::isnan(dens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		//FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				dens[t] = this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			else
			{
				dens[t] = (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			if (std::isnan(dens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_CDFs(double* CDF)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_CDF(this->max_obs+1);
	double dens;

	//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
	// Calculate for j=0
	precomputed_CDF[0] = this->w + (1-this->w) * exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		dens = (1-this->w) * exp( lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp );
		if (std::isnan(dens))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "dens = "<< dens;
			throw nan_detected;
		}
		precomputed_CDF[j] = precomputed_CDF[j-1] + dens;
		if (precomputed_CDF[j] >= 1)
		{
			//FILE_LOG(logDEBUG4) << "CDF >= 1 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
			precomputed_CDF[j] = precomputed_CDF[j-1]; 
		}
	}
	for (int t=0; t<this->T; t++)
	{
		CDF[t] = precomputed_CDF[(int)obs[t]];
		if (std::isnan(CDF[t]))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
			throw nan_detected;
		}
	}
}

void ZiNB::calc_logCDFs(double* logCDF)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_logCDF(this->max_obs+1);
	double logdens;

	//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
	// Calculate for j=0
	precomputed_logCDF[0] = log( this->w + (1-this->w) * exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp ) );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		logdens = log(1-this->w) + lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp;
		if (std::isnan(logdens))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "logdens = "<< logdens;
			throw nan_detected;
		}
		precomputed_logCDF[j] = log( exp(precomputed_logCDF[j-1]) + exp(logdens) );
		if (precomputed_logCDF[j] >= 0)
		{
			//FILE_LOG(logDEBUG4) << "logCDF >= 0 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
			precomputed_logCDF[j] = precomputed_logCDF[j-1]; 
		}
	}
	for (int t=0; t<this->T; t++)
	{
		logCDF[t] = precomputed_logCDF[(int)obs[t]];
		if (std::isnan(logCDF[t]))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
			throw nan_detected;
		}
	}
}

void ZiNB::copy(Density* other)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	ZiNB* o = (ZiNB*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
	this->w = o->w;
}

double ZiNB::getLogDensityAt(int x)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->T; t++)
	{
		mean += obs[t];
	}
	mean = mean / this->T;
	for(int t=0; t<this->T; t++)
	{
		variance += pow(obs[t] - mean, 2);
	}
	variance = variance / this->T;
	// Calculate logdensity
	lGammaR = lgamma(this->size);
	lGammaRplusX = lgamma(this->size + x);
	lxfactorial = this->lxfactorials[x];
	if (x == 0)
	{
		logdens = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp ) );
	}
	else
	{
		logdens = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp;
	}
	if (std::isnan(logdens))
	{
		//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
		//FILE_LOG(logERROR) << "logdens = "<< logdens;
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZiNB::get_mean()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->size*(1-this->prob)/this->prob;
}

double ZiNB::get_variance()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->size*(1-this->prob)/this->prob/this->prob; //TODO: Is this correct?
}

DensityName ZiNB::get_name()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZERO_INFLATED_NEGATIVE_BINOMIAL);
}

double ZiNB::get_size()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->size);
}

double ZiNB::get_prob()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->prob);
}

double ZiNB::get_w()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->w);
}


// ============================================================
// Negative Binomial density
// ============================================================

// Constructor and Destructor ---------------------------------
NegativeBinomial::NegativeBinomial()
{
	this->lxfactorials = NULL;
}

NegativeBinomial::NegativeBinomial(int* observations, int T, double size, double prob)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->size = size;
	this->prob = prob;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	if (this->obs != NULL)
	{
		this->max_obs = intMax(observations, T);
		this->lxfactorials = (double*) Calloc(max_obs+1, double);
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_obs; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

NegativeBinomial::~NegativeBinomial()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->lxfactorials != NULL)
	{
		Free(this->lxfactorials);
	}
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_logdensities(double* logdens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		std::vector<double> logdens_per_read(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			logdens_per_read[j] = lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logdens[t] = logdens_per_read[(int) this->obs[t]];
			//FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (std::isnan(logdens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		//FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			logdens[t] = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			//FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (std::isnan(logdens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_densities(double* dens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		std::vector<double> dens_per_read(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			dens_per_read[j] = exp( lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp );
		}
		for (int t=0; t<this->T; t++)
		{
			dens[t] = dens_per_read[(int) this->obs[t]];
			//FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (std::isnan(dens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		//FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			//FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (std::isnan(dens[t]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_CDFs(double* CDF)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_CDF(this->max_obs+1);
	double dens;

	//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
	// Calculate for j=0
	precomputed_CDF[0] = exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		dens = exp( lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp );
		if (std::isnan(dens))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "dens = "<< dens;
			throw nan_detected;
		}
		precomputed_CDF[j] = precomputed_CDF[j-1] + dens;
		if (precomputed_CDF[j] >= 1)
		{
			//FILE_LOG(logDEBUG4) << "CDF >= 1 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
			precomputed_CDF[j] = precomputed_CDF[j-1]; 
		}
	}
	for (int t=0; t<this->T; t++)
	{
		CDF[t] = precomputed_CDF[(int)obs[t]];
		if (std::isnan(CDF[t]))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
			throw nan_detected;
		}
	}
}

void NegativeBinomial::calc_logCDFs(double* logCDF)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_logCDF(this->max_obs+1);
	double logdens;

	//FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
	// Calculate for j=0
	precomputed_logCDF[0] = lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp;
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		logdens = lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp;
		if (std::isnan(logdens))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "logdens = "<< logdens;
			throw nan_detected;
		}
		precomputed_logCDF[j] = log( exp(precomputed_logCDF[j-1]) + exp(logdens) );
		if (precomputed_logCDF[j] >= 0)
		{
			//FILE_LOG(logDEBUG4) << "logCDF >= 0 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
			precomputed_logCDF[j] = precomputed_logCDF[j-1]; 
		}
	}
	for (int t=0; t<this->T; t++)
	{
		logCDF[t] = precomputed_logCDF[(int)obs[t]];
		if (std::isnan(logCDF[t]))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
			throw nan_detected;
		}
	}
}

void NegativeBinomial::update(double* weights)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	//FILE_LOG(logDEBUG1) << "size = "<<this->size << ", prob = "<<this->prob;
	double eps = 1e-4;
	double kmax = 20;
	double numerator, denominator, size0, DigammaSize, TrigammaSize;
	double F, dFdSize, FdivM;
	double logp = log(this->prob);
	// Update prob (p)
	numerator=denominator=0.0;
// 	clock_t time, dtime;
// 	time = clock();
	for (int t=0; t<this->T; t++)
	{
		numerator += weights[t] * this->size;
		denominator += weights[t] * (this->size + this->obs[t]);
	}
	this->prob = numerator/denominator; // Update this->prob
// 	logp = log(this->prob); // Update of size is done with new prob
	
// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG1) << "updateP(): "<<dtime<< " clicks";

	// Update of size with Newton Method
	size0 = this->size;
// 	time = clock();
	// Select strategy for computing digammas
	if (this->max_obs <= this->T)
	{
		//FILE_LOG(logDEBUG2) << "Precomputing digammas in " << __func__ << " for every obs[t], because max(O)<=T";
		std::vector<double> DigammaSizePlusX(this->max_obs+1);
		std::vector<double> TrigammaSizePlusX(this->max_obs+1);
		for (int k=0; k<kmax; k++)
		{
			F=dFdSize=0.0;
			DigammaSize = digamma(size0); // boost::math::digamma<>(size0);
			TrigammaSize = trigamma(size0); // boost::math::digamma<>(size0);
			// Precompute the digammas by iterating over all possible values of the observation vector
			for (int j=0; j<=this->max_obs; j++)
			{
				DigammaSizePlusX[j] = digamma(size0+j);
				TrigammaSizePlusX[j] = trigamma(size0+j);
			}
			for(int t=0; t<this->T; t++)
			{
				if(this->obs[t]==0)
				{
					F += weights[t] * logp;
					//dFdSize+=0;
				}
				if(this->obs[t]!=0)
				{
					F += weights[t] * (logp - DigammaSize + DigammaSizePlusX[(int)obs[t]]);
					dFdSize += weights[t] * (-TrigammaSize + TrigammaSizePlusX[(int)obs[t]]);
				}
			}
			FdivM = F/dFdSize;
// Rprintf("k = %d, F = %g, dFdSize = %g, FdivM = %g, size0 = %g\n", k, F, dFdSize, FdivM, size0);
			if (FdivM < size0)
			{
				size0 = size0-FdivM;
			}
			else if (FdivM >= size0)
			{
				size0 = size0/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
		}
	}
	else
	{
		//FILE_LOG(logDEBUG2) << "Computing digammas in " << __func__ << " for every t, because max(O)>T";
		double DigammaSizePlusX, TrigammaSizePlusX;
		for (int k=0; k<kmax; k++)
		{
			F = dFdSize = 0.0;
			DigammaSize = digamma(size0); // boost::math::digamma<>(size0);
			TrigammaSize = trigamma(size0); // boost::math::digamma<>(size0);
			for(int t=0; t<this->T; t++)
			{
				DigammaSizePlusX = digamma(size0+this->obs[t]); //boost::math::digamma<>(size0+this->obs[ti]);
				TrigammaSizePlusX = trigamma(size0+this->obs[t]);
				if(this->obs[t]==0)
				{
					F += weights[t] * logp;
					//dFdSize+=0;
				}
				if(this->obs[t]!=0)
				{
					F += weights[t] * (logp - DigammaSize + DigammaSizePlusX);
					dFdSize += weights[t] * (-TrigammaSize + TrigammaSizePlusX);
				}
			}
			FdivM = F/dFdSize;
			if (FdivM < size0)
			{
				size0 = size0-FdivM;
			}
			else if (FdivM >= size0)
			{
				size0 = size0/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
		}
	}
	this->size = size0;
	//FILE_LOG(logDEBUG1) << "size = "<<this->size << ", prob = "<<this->prob;

// 	dtime = clock() - time;
// 	//FILE_LOG(logDEBUG1) << "updateR(): "<<dtime<< " clicks";

}

// void NegativeBinomial::update(double* weight)
// {
// 	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	double eps = 1e-4, kmax;
// 	double numerator, denominator, rhere, dr, Fr, dFrdr, DigammaR, DigammaRplusDR;
// 	// Update p
// 	numerator=denominator=0.0;
// // 	clock_t time, dtime;
// // 	time = clock();
// 	for (int t=0; t<this->T; t++)
// 	{
// 		numerator+=weight[t]*this->size;
// 		denominator+=weight[t]*(this->size+this->obs[t]);
// 	}
// 	this->prob = numerator/denominator; // Update of r is now done with updated p
// 	double logp = log(this->prob);
// // 	dtime = clock() - time;
// // 	//FILE_LOG(logDEBUG1) << "updateP(): "<<dtime<< " clicks";
// 	// Update of r with Newton Method
// 	rhere = this->size;
// 	dr = 0.00001;
// 	kmax = 20;
// // 	time = clock();
// 	// Select strategy for computing digammas
// 	if (this->max_obs <= this->T)
// 	{
// 		//FILE_LOG(logDEBUG3) << "Precomputing digammas in " << __func__ << " for every obs[t], because max(O)<=T";
// 		std::vector<double> DigammaRplusX(this->max_obs+1);
// 		std::vector<double> DigammaRplusDRplusX(this->max_obs+1);
// 		for (int k=1; k<kmax; k++)
// 		{
// 			Fr=dFrdr=0.0;
// 			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
// 			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
// 			// Precompute the digammas by iterating over all possible values of the observation vector
// 			for (int j=0; j<=this->max_obs; j++)
// 			{
// 				DigammaRplusX[j] = digamma(rhere+j);
// 				DigammaRplusDRplusX[j] = digamma(rhere+dr+j);
// 			}
// 			for(int t=0; t<this->T; t++)
// 			{
// 				if(this->obs[t]==0)
// 				{
// 					Fr+=weight[t]*logp;
// 					//dFrdr+=0;
// 				}
// 				if(this->obs[t]!=0)
// 				{
// 					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX[(int)obs[t]]);
// 					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX[(int)obs[t]]-DigammaRplusX[(int)obs[t]]);
// 				}
// 			}
// 			if(fabs(Fr)<eps)
// 			{
// 				break;
// 			}
// 			if(Fr/dFrdr<rhere)
// 			{
// 				rhere=rhere-Fr/dFrdr;
// 			}
// 			else if (Fr/dFrdr>=rhere)
// 			{
// 				rhere=rhere/2.0;
// 			}
// 		}
// 	}
// 	else
// 	{
// 		//FILE_LOG(logDEBUG2) << "Computing digammas in " << __func__ << " for every t, because max(O)>T";
// 		double DigammaRplusX, DigammaRplusDRplusX;
// 		for (int k=1; k<kmax; k++)
// 		{
// 			Fr = dFrdr = 0.0;
// 			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
// 			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
// 			for(int t=0; t<this->T; t++)
// 			{
// 				DigammaRplusX = digamma(rhere+this->obs[t]); //boost::math::digamma<>(rhere+this->obs[ti]);
// 				DigammaRplusDRplusX = digamma(rhere+dr+this->obs[t]); // boost::math::digamma<>(rhere+dr+this->obs[ti]);
// 				if(this->obs[t]==0)
// 				{
// 					Fr+=weight[t]*logp;
// 					//dFrdr+=0;
// 				}
// 				if(this->obs[t]!=0)
// 				{
// 					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX);
// 					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
// 				}
// 			}
// 			if(fabs(Fr)<eps)
// 			{
// 				break;
// 			}
// 			if(Fr/dFrdr<rhere)
// 			{
// 				rhere=rhere-Fr/dFrdr;
// 			}
// 			else if (Fr/dFrdr>=rhere)
// 			{
// 				rhere=rhere/2.0;
// 			}
// 		}
// 	}
// 	this->size = rhere;
// 	//FILE_LOG(logDEBUG1) << "r = "<<this->size << ", p = "<<this->prob;
// 
// // 	dtime = clock() - time;
// // 	//FILE_LOG(logDEBUG1) << "updateR(): "<<dtime<< " clicks";
// 
// }

void NegativeBinomial::copy(Density* other)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	NegativeBinomial* o = (NegativeBinomial*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
}

double NegativeBinomial::getLogDensityAt(int x)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->T; t++)
	{
		mean += obs[t];
	}
	mean = mean / this->T;
	for(int t=0; t<this->T; t++)
	{
		variance += pow(obs[t] - mean, 2);
	}
	variance = variance / this->T;
	// Calculate logdensity
	lGammaR = lgamma(this->size);
	lGammaRplusX = lgamma(this->size + x);
	lxfactorial = this->lxfactorials[x];
	logdens = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp;
	if (std::isnan(logdens))
	{
		//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
		//FILE_LOG(logERROR) << "logdens = "<< logdens;
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double NegativeBinomial::get_mean()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob;
}

double NegativeBinomial::get_variance()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob/this->prob;
}

DensityName NegativeBinomial::get_name()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NEGATIVE_BINOMIAL);
}

double NegativeBinomial::get_size()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->size);
}

double NegativeBinomial::get_prob()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->prob);
}


// ============================================================
// Zero Inflation density
// ============================================================

// Constructor and Destructor ---------------------------------
ZeroInflation::ZeroInflation() {}

ZeroInflation::ZeroInflation(int* observations, int T)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
}

ZeroInflation::~ZeroInflation()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void ZeroInflation::calc_logdensities(double* logdens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(obs[t]==0)
		{
			logdens[t] = 0.0;
		};
		if(obs[t]>0)
		{
			logdens[t] = -INFINITY;
		}
		//FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
	}
}

void ZeroInflation::calc_densities(double* dens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(obs[t]==0)
		{
			dens[t] = 1.0;
		}
		if(obs[t]>0)
		{
			dens[t] = 0.0;
		}
		//FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
	}
}

void ZeroInflation::copy(Density*) {}

void ZeroInflation::update(double*)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

double ZeroInflation::getLogDensityAt(int x)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logdens;
	// Calculate logdensity
	if (x == 0)
	{
		logdens = 0;
	}
	else
	{
		logdens = -INFINITY;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZeroInflation::get_mean()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

double ZeroInflation::get_variance()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

DensityName ZeroInflation::get_name()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZERO_INFLATION);
}

	
// ============================================================
// Multivariate Copula Approximation
// ============================================================

// Constructor and Destructor ---------------------------------
MVCopulaApproximation::MVCopulaApproximation(int** multiobservations, int T, std::vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->multi_obs = multiobservations;
	this->T = T;
	// these are the marginal distributions (we need their CDF function)
	this->marginals = marginals;
	this->Nmod = this->marginals.size();
	this->cor_matrix_inv = cor_matrix_inv;
	this->cor_matrix_determinant = cor_matrix_determinant;
}

MVCopulaApproximation::~MVCopulaApproximation()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int imod=0; imod<this->Nmod; imod++)
	{
		delete this->marginals[imod];
	}
}

// Methods ----------------------------------------------------
void MVCopulaApproximation::calc_logdensities(double* logdens)
{
// Rprintf("new state\n");
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	// Calculate logdensities for marginals
	double** marginals_logdensities = CallocDoubleMatrix(this->Nmod, this->T);
	double** marginals_CDFs = CallocDoubleMatrix(this->Nmod, this->T);
	for (int imod=0; imod<this->Nmod; imod++)
	{
		//FILE_LOG(logDEBUG2) << __func__ << ": calculating marginals for imod = " << imod;
		this->marginals[imod]->calc_logdensities(marginals_logdensities[imod]);
		this->marginals[imod]->calc_CDFs(marginals_CDFs[imod]);
	}
	// Calculate multivariate Copula approximation
	//FILE_LOG(logDEBUG2) << __func__ << ": calculate Copula approximation";
	double sum, uniform, exponent, exponentTemp;
	double* z = (double*) Calloc(this->Nmod, double);
	for (int t=0; t<this->T; t++)
	{
		sum = 0.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			sum += marginals_logdensities[imod][t];
			uniform = marginals_CDFs[imod][t];
			z[imod] = qnorm(uniform, 0, 1, 1, 0);
			if (std::isnan(z[imod]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "uniform = "<< uniform;
				//FILE_LOG(logERROR) << "z[imod] = "<< z[imod];
				throw nan_detected;
			}
// if (t==0)
// {
// 	Rprintf("\nmarginal_logdensities[imod=%d][%d] = %g\n", imod, t, marginals_logdensities[imod][t]);
// 	Rprintf("sum = %g\n", sum);
// 	Rprintf("uniform = %g\n", uniform);
// 	Rprintf("z[imod=%d] = %g\n", imod, z[imod]);
// }
		}
		exponent = 0.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			exponentTemp = 0.0;
			for(int jmod=0; jmod<Nmod; jmod++)
			{
				if(imod==jmod)
				{
					exponentTemp += z[jmod] * (this->cor_matrix_inv[imod * Nmod + jmod] - 1);
				}
				else
				{
					exponentTemp += z[jmod] * this->cor_matrix_inv[imod * Nmod + jmod];
				}
				if (std::isnan(exponentTemp))
				{
					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					//FILE_LOG(logERROR) << "exponentTemp = "<< exponentTemp;
					//FILE_LOG(logERROR) << "cor_matrix_inv = "<< cor_matrix_inv[imod * Nmod + jmod];
					//FILE_LOG(logERROR) << "z["<<jmod<<"] = "<< z[jmod];
					throw nan_detected;
				}
			}
			exponent += exponentTemp * z[imod];
			if (std::isnan(exponent))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				//FILE_LOG(logERROR) << "exponentTemp = "<< exponentTemp;
				//FILE_LOG(logERROR) << "z["<<imod<<"] = "<< z[imod];
				//FILE_LOG(logERROR) << "exponent = "<< exponent;
				throw nan_detected;
			}
		}
// 		logdens[t] = log(1/sqrt(this->cor_matrix_determinant)) - 0.5 * exponent + sum;
		logdens[t] = -0.5 * log(this->cor_matrix_determinant) - 0.5 * exponent + sum;
		if (std::isnan(logdens[t]))
		{
			//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			//FILE_LOG(logERROR) << "cor_matrix_determinant = " << this->cor_matrix_determinant;
			//FILE_LOG(logERROR) << "sum = " << sum;
			//FILE_LOG(logERROR) << "exponentTemp = " << exponentTemp;
			//FILE_LOG(logERROR) << "exponent = " << exponent;
			//FILE_LOG(logERROR) << "logdens["<<t<<"] = " << logdens[t];
			throw nan_detected;
		}		
// if (t==0)
// {
// 	Rprintf("\nlogdens[%d] = %g\n", t, logdens[t]);
// 	Rprintf("-0.5*exponent = %g\n", -0.5*exponent);
// 	Rprintf("sum = %g\n", sum);
// 	Rprintf("cor_matrix_determinant = %g\n", this->cor_matrix_determinant);
// 	Rprintf("-0.5*log(cor_matrix_determinant) = %g\n", -0.5*log(this->cor_matrix_determinant));
// }
	}

	// Clean up
	FreeDoubleMatrix(marginals_logdensities, this->Nmod);
	FreeDoubleMatrix(marginals_CDFs, this->Nmod);
	Free(z);
}

void MVCopulaApproximation::calc_densities(double* dens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->calc_logdensities(dens);

	for (int t=0; t<this->T; t++)
	{
		dens[t] = exp( dens[t] );
	}
}

// Getter and Setter ------------------------------------------
DensityName MVCopulaApproximation::get_name()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(OTHER);
}

	
// ============================================================
// Multivariate Product of Bernoullis
// ============================================================

// Constructor and Destructor ---------------------------------
BernoulliProduct::BernoulliProduct(double** multiobservations, bool* binary_states, int T, int Nmod)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->multi_obs = multiobservations;
	this->binary_states = binary_states;
	this->T = T;
	this->Nmod = Nmod;
}

BernoulliProduct::~BernoulliProduct()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void BernoulliProduct::calc_logdensities(double* logdens)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double d, mult;
	double** tempPost = CallocDoubleMatrix(this->Nmod, this->T);

	for (int t=0; t<this->T; t++)
	{
		d = 1.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			//if state[iN] is such that modification imod is unmodified, multiProb[t][imod] is the univariate posterior of being unmodified. 
			//if state[iN] is such that modification imod is modified, multiProb[t][imod] is the univariate posterior of being modified
			if (binary_states[imod])
			{
				mult = 1-this->multi_obs[imod][t];
			}
			else
			{
				mult = this->multi_obs[imod][t];
			}
			if(mult>=1) mult=0.9999999999999;
			if(mult<=0) mult=0.0000000000001;
			d=d*mult;
		}
		logdens[t] = log(d);
	}
	FreeDoubleMatrix(tempPost, this->Nmod);
}

// Getter and Setter ------------------------------------------
DensityName BernoulliProduct::get_name()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(OTHER);
}

