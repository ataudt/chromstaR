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
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->prob = prob;
	this->size = size;
	this->w = w;
	this->lxfactorials = NULL;
	if (this->obs != NULL)
	{
		this->max_obs = intMax(observations, T);
		this->lxfactorials = (double*) calloc(max_obs+1, sizeof(double));
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_obs; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

ZiNB::~ZiNB()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->obs != NULL)
	{
		free(this->lxfactorials);
	}
}

// Methods ----------------------------------------------------
void ZiNB::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGammaRplusX[this->max_obs+1];
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
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
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
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGammaRplusX[this->max_obs+1];
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
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
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
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_CDFs(double* CDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->prob);
	double lGammaR = lgamma(this->size);
	double lppowerr = this->size * log(this->prob);
	// No selection strategy here, because we must precompute CDF to deal with 1s by shifting them
// 	// Select strategy for computing gammas
// 	if (this->max_obs <= this->T)
//	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_obs+1], lGamma2plusX[this->max_obs+1], lHyper[this->max_obs+1], lppowert[this->max_obs+1];
		double precomputed_CDF[this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->size + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->size + j, 2 + j, 1-this->prob));
			lppowert[j] = (1+j) * log1minusp;
			precomputed_CDF[j] = 1 - exp( log(1-this->w) + lppowerr + lppowert[j] + lGamma1plusRplusX[j] + lHyper[j] - lGammaR - lGamma2plusX[j] );
			if (precomputed_CDF[j] == 1)
			{
				FILE_LOG(logDEBUG4) << "CDF = 1 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
				precomputed_CDF[j] = precomputed_CDF[j-1]; 
			}
		}
		for (int t=0; t<this->T; t++)
		{
			CDF[t] = precomputed_CDF[(int)obs[t]];
			if (isnan(CDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
				throw nan_detected;
			}
		}
// 	}
// 	else
// 	{
// 		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
// 		double lGamma1plusRplusX, lGamma2plusX, lHyper, lppowert;
// 		for (int t=0; t<this->T; t++)
//		{
// 			lGamma1plusRplusX = lgamma(1 + this->size + this->obs[t]);
// 			lGamma2plusX = lgamma(2 + this->obs[t]);
// 			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->size + this->obs[t], 2 + this->obs[t], 1-this->prob));
// 			lppowert = (1+this->obs[t]) * log1minusp;
// 			CDF[t] = 1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX ); //TODO: Check formula for log
// 			if(CDF[t] == 0)
// 			{
// 				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t]; //TODO: Check if this works
// // 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
// 			}
// 			if (isnan(CDF[t]))
// 			{
// 				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
// 				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
// 				throw nan_detected;
// 			}
// 		}
// 	}
}

void ZiNB::calc_logCDFs(double* logCDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->prob);
	double lGamma1plusRplusX, lHyper, lGammaR, lGamma2plusX, lppowert, lppowerr;
	lGammaR = lgamma(this->size);
	lppowerr = this->size * log(this->prob);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_obs+1], lGamma2plusX[this->max_obs+1], lHyper[this->max_obs+1], lppowert[this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->size + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->size + j, 2 + j, 1-this->prob));
			lppowert[j] = (1+j) * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logCDF[t] = log(1 - exp( log(1-this->w) + lppowerr + lppowert[(int)obs[t]] + lGamma1plusRplusX[(int)obs[t]] + lHyper[(int)obs[t]] - lGammaR - lGamma2plusX[(int)obs[t]] ));
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0";
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGamma1plusRplusX = lgamma(1 + this->size + this->obs[t]);
			lGamma2plusX = lgamma(2 + this->obs[t]);
			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->size + this->obs[t], 2 + this->obs[t], 1-this->prob));
			lppowert = (1+this->obs[t]) * log1minusp;
			logCDF[t] = log(1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )); //TODO: Check formula for log
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0"; //TODO: Check if this works
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
				throw nan_detected;
			}
		}
	}
}

void ZiNB::copy(Density* other)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	ZiNB* o = (ZiNB*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
	this->w = o->w;
}

double ZiNB::getLogDensityAt(int x)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
	if (isnan(logdens))
	{
		FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
		FILE_LOG(logERROR) << "logdens = "<< logdens;
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZiNB::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->size*(1-this->prob)/this->prob;
}

double ZiNB::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->size*(1-this->prob)/this->prob/this->prob; //TODO: Is this correct?
}

DensityName ZiNB::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZERO_INFLATED_NEGATIVE_BINOMIAL);
}

double ZiNB::get_size()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->size);
}

double ZiNB::get_prob()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->prob);
}

double ZiNB::get_w()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->size = size;
	this->prob = prob;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	if (this->obs != NULL)
	{
		this->max_obs = intMax(observations, T);
		this->lxfactorials = (double*) calloc(max_obs+1, sizeof(double));
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_obs; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

NegativeBinomial::~NegativeBinomial()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->lxfactorials != NULL)
	{
		free(this->lxfactorials);
	}
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double logdens_per_read [this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			logdens_per_read[j] = lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logdens[t] = logdens_per_read[(int) this->obs[t]];
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			logdens[t] = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double dens_per_read [this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			dens_per_read[j] = exp( lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp );
		}
		for (int t=0; t<this->T; t++)
		{
			dens[t] = dens_per_read[(int) this->obs[t]];
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_CDFs(double* CDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->prob);
	double lGammaR = lgamma(this->size);
	double lppowerr = this->size * log(this->prob);
	// No selection strategy here, because we must precompute CDF to deal with 1s by shifting them
// 	// Select strategy for computing gammas
// 	if (this->max_obs <= this->T)
// 	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_obs+1], lGamma2plusX[this->max_obs+1], lHyper[this->max_obs+1], lppowert[this->max_obs+1];
		double precomputed_CDF[this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->size + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->size + j, 2 + j, 1-this->prob));
			lppowert[j] = (1+j) * log1minusp;
			precomputed_CDF[j] = 1 - exp( lppowerr + lppowert[j] + lGamma1plusRplusX[j] + lHyper[j] - lGammaR - lGamma2plusX[j] );
			if (precomputed_CDF[j] == 1)
			{
				FILE_LOG(logDEBUG4) << "CDF = 1 for obs[t] = "<<j<< ", shifting to value of obs[t] = "<<j-1;
				precomputed_CDF[j] = precomputed_CDF[j-1]; 
			}
		}
		for (int t=0; t<this->T; t++)
		{
			CDF[t] = precomputed_CDF[(int)obs[t]];
			if (isnan(CDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
				throw nan_detected;
			}
		}
// 	}
// 	else
// 	{
// 		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
// 		double lGamma1plusRplusX, lGamma2plusX, lHyper, lppowert;
// 		for (int t=0; t<this->T; t++)
// 		{
// 			lGamma1plusRplusX = lgamma(1 + this->size + this->obs[t]);
// 			lGamma2plusX = lgamma(2 + this->obs[t]);
// 			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->size + this->obs[t], 2 + this->obs[t], 1-this->prob));
// 			lppowert = (1+this->obs[t]) * log1minusp;
// 			CDF[t] = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX ); //TODO: Check formula for log
// 			if(CDF[t] == 1)
// 			{
// 				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t]; //TODO: Check if this works
// // 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
// 			}
// 			if (isnan(CDF[t]))
// 			{
// 				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
// 				FILE_LOG(logERROR) << "CDF["<<t<<"] = "<< CDF[t];
// 				throw nan_detected;
// 			}
// 		}
// 	}
}

void NegativeBinomial::calc_logCDFs(double* logCDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->prob);
	double lGamma1plusRplusX, lHyper, lGammaR, lGamma2plusX, lppowert, lppowerr;
	lGammaR = lgamma(this->size);
	lppowerr = this->size * log(this->prob);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_obs+1], lGamma2plusX[this->max_obs+1], lHyper[this->max_obs+1], lppowert[this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->size + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->size + j, 2 + j, 1-this->prob));
			lppowert[j] = (1+j) * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logCDF[t] = log(1 - exp( lppowerr + lppowert[(int)obs[t]] + lGamma1plusRplusX[(int)obs[t]] + lHyper[(int)obs[t]] - lGammaR - lGamma2plusX[(int)obs[t]] ));
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0";
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGamma1plusRplusX = lgamma(1 + this->size + this->obs[t]);
			lGamma2plusX = lgamma(2 + this->obs[t]);
			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->size + this->obs[t], 2 + this->obs[t], 1-this->prob));
			lppowert = (1+this->obs[t]) * log1minusp;
			logCDF[t] = log(1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )); //TODO: Check formula for log
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 1"; //TODO: Check if this works
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logCDF["<<t<<"] = "<< logCDF[t];
				throw nan_detected;
			}
		}
	}
}

void NegativeBinomial::update(double* weight)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double eps = 1e-4, kmax;
	double numerator, denominator, rhere, dr, Fr, dFrdr, DigammaR, DigammaRplusDR;
	// Update p
	numerator=denominator=0.0;
// 	clock_t time, dtime;
// 	time = clock();
	for (int t=0; t<this->T; t++)
	{
		numerator+=weight[t]*this->size;
		denominator+=weight[t]*(this->size+this->obs[t]);
	}
	this->prob = numerator/denominator; // Update of r is now done with updated p
	double logp = log(this->prob);
// 	dtime = clock() - time;
// 	FILE_LOG(logDEBUG1) << "updateP(): "<<dtime<< " clicks";
	// Update of r with Newton Method
	rhere = this->size;
	dr = 0.00001;
	kmax = 20;
// 	time = clock();
	// Select strategy for computing digammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG3) << "Precomputing digammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double DigammaRplusX[this->max_obs+1], DigammaRplusDRplusX[this->max_obs+1];
		for (int k=1; k<kmax; k++)
		{
			Fr=dFrdr=0.0;
			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
			// Precompute the digammas by iterating over all possible values of the observation vector
			for (int j=0; j<=this->max_obs; j++)
			{
				DigammaRplusX[j] = digamma(rhere+j);
				DigammaRplusDRplusX[j] = digamma(rhere+dr+j);
			}
			for(int t=0; t<this->T; t++)
			{
				if(this->obs[t]==0)
				{
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX[(int)obs[t]]);
					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX[(int)obs[t]]-DigammaRplusX[(int)obs[t]]);
				}
			}
			if(fabs(Fr)<eps)
{
				break;
			}
			if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
			if(Fr/dFrdr>rhere) rhere=rhere/2.0;
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing digammas in " << __func__ << " for every t, because max(O)>T";
		double DigammaRplusX, DigammaRplusDRplusX;
		for (int k=1; k<kmax; k++)
		{
			Fr = dFrdr = 0.0;
			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
			for(int t=0; t<this->T; t++)
			{
				DigammaRplusX = digamma(rhere+this->obs[t]); //boost::math::digamma<>(rhere+this->obs[ti]);
				DigammaRplusDRplusX = digamma(rhere+dr+this->obs[t]); // boost::math::digamma<>(rhere+dr+this->obs[ti]);
				if(this->obs[t]==0)
				{
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX);
					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
				}
			}
			if(fabs(Fr)<eps)
			{
				break;
			}
			if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
			if(Fr/dFrdr>rhere) rhere=rhere/2.0;
		}
	}
	this->size = rhere;
	FILE_LOG(logDEBUG1) << "r = "<<this->size << ", p = "<<this->prob;

// 	dtime = clock() - time;
// 	FILE_LOG(logDEBUG1) << "updateR(): "<<dtime<< " clicks";

}

void NegativeBinomial::copy(Density* other)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	NegativeBinomial* o = (NegativeBinomial*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
}

double NegativeBinomial::getLogDensityAt(int x)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
	if (isnan(logdens))
	{
		FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
		FILE_LOG(logERROR) << "logdens = "<< logdens;
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double NegativeBinomial::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob;
}

double NegativeBinomial::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob/this->prob;
}

DensityName NegativeBinomial::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NEGATIVE_BINOMIAL);
}

double NegativeBinomial::get_size()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->size);
}

double NegativeBinomial::get_prob()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->prob);
}


// ============================================================
// Zero Inflation density
// ============================================================

// Constructor and Destructor ---------------------------------
ZeroInflation::ZeroInflation() {}

ZeroInflation::ZeroInflation(int* observations, int T)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
}

ZeroInflation::~ZeroInflation()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void ZeroInflation::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
		FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
	}
}

void ZeroInflation::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
		FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
	}
}

void ZeroInflation::copy(Density* other) {}

void ZeroInflation::update(double* weight)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

double ZeroInflation::getLogDensityAt(int x)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
	if (x == 0)
	{
		logdens = 0;
	}
	else
	{
		logdens = -100;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZeroInflation::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

double ZeroInflation::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

DensityName ZeroInflation::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZERO_INFLATION);
}

	
// ============================================================
// Multivariate Copula Approximation
// ============================================================

// Constructor and Destructor ---------------------------------
MVCopulaApproximation::MVCopulaApproximation(int** multiobservations, int T, std::vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
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
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int imod; imod<this->Nmod; imod++)
	{
		delete this->marginals[imod];
	}
}

// Methods ----------------------------------------------------
void MVCopulaApproximation::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	// Calculate logdensities for marginals
	double** marginals_logdensities = allocDoubleMatrix(this->Nmod, this->T);
	double** marginals_CDFs = allocDoubleMatrix(this->Nmod, this->T);
	for (int imod=0; imod<this->Nmod; imod++)
	{
		FILE_LOG(logDEBUG2) << __func__ << ": calculating marginals for imod = " << imod;
		this->marginals[imod]->calc_logdensities(marginals_logdensities[imod]);
		this->marginals[imod]->calc_CDFs(marginals_CDFs[imod]);
	}
	// Calculate multivariate Copula approximation
	FILE_LOG(logDEBUG2) << __func__ << ": calculate Copula approximation";
	double sum, uniform, exponent, exponentTemp;
	double* z = (double*) calloc(this->Nmod, sizeof(double));
	for (int t=0; t<this->T; t++)
	{
		sum = 0.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			sum += marginals_logdensities[imod][t];
			uniform = marginals_CDFs[imod][t];
			z[imod] = qnorm(uniform, 0, 1, 1, 0);
			if (isnan(z[imod]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "uniform = "<< uniform;
				FILE_LOG(logERROR) << "z[imod] = "<< z[imod];
				throw nan_detected;
			}
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
				if (isnan(exponentTemp))
				{
					FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					FILE_LOG(logERROR) << "exponentTemp = "<< exponentTemp;
					FILE_LOG(logERROR) << "cor_matrix_inv = "<< cor_matrix_inv[imod * Nmod + jmod];
					FILE_LOG(logERROR) << "z["<<jmod<<"] = "<< z[jmod];
					throw nan_detected;
				}
			}
			exponent += exponentTemp * z[imod];
			if (isnan(exponent))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "exponentTemp = "<< exponentTemp;
				FILE_LOG(logERROR) << "z["<<imod<<"] = "<< z[imod];
				FILE_LOG(logERROR) << "exponent = "<< exponent;
				throw nan_detected;
			}
		}
// 		logdens[t] = log(1/sqrt(this->cor_matrix_determinant)) - 0.5 * exponent + sum;
		logdens[t] = -0.5 * log(this->cor_matrix_determinant) - 0.5 * exponent + sum;
		if (isnan(logdens[t]))
		{
			FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
			FILE_LOG(logERROR) << "cor_matrix_determinant = " << this->cor_matrix_determinant;
			FILE_LOG(logERROR) << "sum = " << sum;
			FILE_LOG(logERROR) << "exponentTemp = " << exponentTemp;
			FILE_LOG(logERROR) << "exponent = " << exponent;
			FILE_LOG(logERROR) << "logdens["<<t<<"] = " << logdens[t];
			throw nan_detected;
		}		
	}

	// Clean up
	freeDoubleMatrix(marginals_logdensities, this->Nmod);
	freeDoubleMatrix(marginals_CDFs, this->Nmod);
	free(z);
}

void MVCopulaApproximation::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->calc_logdensities(dens);

	for (int t=0; t<this->T; t++)
	{
		dens[t] = exp( dens[t] );
	}
}

// Getter and Setter ------------------------------------------
DensityName MVCopulaApproximation::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(OTHER);
}

	
// ============================================================
// Multivariate Product of Bernoullis
// ============================================================

// Constructor and Destructor ---------------------------------
BernoulliProduct::BernoulliProduct(double** multiobservations, bool* binary_states, int T, int Nmod)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->multi_obs = multiobservations;
	this->binary_states = binary_states;
	this->T = T;
	this->Nmod = Nmod;
}

BernoulliProduct::~BernoulliProduct()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void BernoulliProduct::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double d, mult;
	double** tempPost = allocDoubleMatrix(this->Nmod, this->T);

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
	freeDoubleMatrix(tempPost, this->Nmod);
}

// Getter and Setter ------------------------------------------
DensityName BernoulliProduct::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(OTHER);
}

