#include <Rmath.h> // digamma() and qnorm()
#include <vector> // storing density functions in MVCopula
#include "utility.h" // FILE_LOG(), intMax(), allocDoubleMatrix()

#ifndef DENSITIES_H
#define DENSITIES_H

enum whichvariate {UNIVARIATE, MULTIVARIATE};
enum DensityName {ZERO_INFLATION, NEGATIVE_BINOMIAL, ZERO_INFLATED_NEGATIVE_BINOMIAL, OTHER};

class Density {
	public:
		// Constructor and Destructor
		virtual ~Density() {};
		// Methods
		virtual void calc_logdensities(double*) {};
		virtual void calc_densities(double*) {};
		virtual void calc_logCDFs(double*) {};
		virtual void calc_CDFs(double*) {};
		virtual void update(double*) {}; 
		virtual void copy(Density*) {};
		virtual double getLogDensityAt(int) { return(0); };
		// Getter and Setter
		virtual DensityName get_name() { return(OTHER); };
		virtual double get_mean() { return(0); };
		virtual double get_variance() { return(0); };

};  


class ZiNB : public Density {
	public:
		// Constructor and Destructor
		ZiNB();
		ZiNB(int* observations, int T, double size, double prob, double w);
		~ZiNB();
	
		// Methods
		void calc_logdensities(double* logdensity);
		void calc_densities(double* logdensity);
		void calc_logCDFs(double* logCDF);
		void calc_CDFs(double* CDF);
		void copy(Density* other);
		double getLogDensityAt(int x);
	
		// Getter and Setter
		double get_mean();
		double get_variance();
		DensityName get_name();
		double get_size();
		double get_prob();
		double get_w();

	private:
		// Member variables
		double size; ///< parameter of the distribution
		double prob; ///< parameter of the distribution
		double w; ///< parameter of the distribution
		int* obs; ///< vector [T] of observations
		int T; ///< length of observation vector
// 		double* weight; ///< temporary storage for weights in update()
		int max_obs; ///< maximum value in *obs
		double* lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};

class NegativeBinomial : public Density {
	public:
		// Constructor and Destructor
		NegativeBinomial();
		NegativeBinomial(int* observations, int T, double size, double prob);
		~NegativeBinomial();

		// Methods
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);
		void calc_logCDFs(double* logCDF);
		void calc_CDFs(double* CDF);
		void update(double* weight);
		void copy(Density* other);
		double getLogDensityAt(int x);

		// Getter and Setter
		double get_mean();
		double get_variance();
		DensityName get_name();
		double get_size();
		double get_prob();

	private:
		// Member variables
		double size; ///< parameter of the distribution
		double prob; ///< parameter of the distribution
		int* obs; ///< vector [T] of observations
		int T; ///< length of observation vector
		int max_obs; ///< maximum value in *obs
		double* lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};


class ZeroInflation : public Density {
	public:
		// Constructor and Destructor
		ZeroInflation();
		ZeroInflation(int* observations, int T);
		~ZeroInflation();

		// Methods
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);
		void update(double* weight);
		void copy(Density* other);
		double getLogDensityAt(int x);

		// Getters and Setters
		double get_mean();
		double get_variance();
		DensityName get_name();

	private:
		// Member variables
		int* obs; ///< vector [T] of observations
		int T; ///< length of observation vector
};


class MVCopulaApproximation : public Density {
	public:
		// Constructor and Destructor
		MVCopulaApproximation(int** multiobservations, int T, std::vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant);
		~MVCopulaApproximation();
	
		// Methods
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);

		// Getters and Setters
		DensityName get_name();

	private:
		// Member variables
		int Nmod; ///< number of modifications
		int** multi_obs; ///< matrix [Nmod x T] of observations
		int T; ///< length of observation vector
		std::vector<Density*> marginals; ///< vector [Nmod] of marginal distributions
		double* cor_matrix_inv; ///< vector with elements of the inverse of the correlation matrix
		double cor_matrix_determinant; ///< determinant of the correlation matrix
};


class BernoulliProduct : public Density {
	public:
		// Constructor and Destructor
		BernoulliProduct(double** multiobservations, bool* binary_states, int T, int Nmod);
		~BernoulliProduct();
		// Methods
		void calc_logdensities(double* logdens);
		// Getters and Setters
		DensityName get_name();

	private:
		// Member variables
		double** multi_obs;
		bool* binary_states;
		int T;
		int Nmod;
};


#endif
