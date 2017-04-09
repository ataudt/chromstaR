


#include "utility.h"
#include "scalehmm.h"
#include <vector> // storing density functions in multivariate
#include <string> // strcmp
#include <Rmath.h> //runif

#ifdef _OPENMP
#include <omp.h>
#endif
// #if defined TARGET_OS_MAC || defined __APPLE__
// #include <libiomp/omp.h> // parallelization options on mac
// #elif defined __linux__ || defined _WIN32 || defined _WIN64
// #include <omp.h> // parallelization options
// #endif

extern "C"
void univariate_hmm(int* O, int* T, int* N, double* size, double* prob, int* maxiter, int* maxtime, double* eps, double* posteriors, double* densities, bool* keep_densities, double* A, double* proba, double* loglik, double* weights, int* iniproc, double* initial_size, double* initial_prob, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff, int* verbosity);

extern "C"
void multivariate_hmm(int* O, int* T, int* N, int *Nmod, double* comb_states, double* size, double* prob, double* w, double* cor_matrix_inv, double* det, int* maxiter, int* maxtime, double* eps, double* posteriors, bool* keep_posteriors, double* densities, bool* keep_densities, int* states, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* verbosity);

extern "C"
void univariate_cleanup();

extern "C"
void multivariate_cleanup(int* Nmod);

extern "C"
void array3D_which_max(double* array3D, int* dim, int* ind_max);

extern "C"
void array2D_mean(double* array2D, int* dim, double* mean);

extern "C"
void array3D_mean(double* array3D, int* dim, double* mean);