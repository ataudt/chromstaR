#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "R_interface.h"


R_NativePrimitiveArgType arg1[] = {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, LGLSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP};
R_NativePrimitiveArgType arg2[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, LGLSXP, REALSXP, LGLSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP};
R_NativePrimitiveArgType arg4[] = {INTSXP};
R_NativePrimitiveArgType arg5[] = {REALSXP, INTSXP, INTSXP};
R_NativePrimitiveArgType arg6[] = {REALSXP, INTSXP, REALSXP};
R_NativePrimitiveArgType arg7[] = {REALSXP, INTSXP, REALSXP};

static const R_CMethodDef CEntries[]  = {
    {"C_univariate_hmm", (DL_FUNC) &univariate_hmm, 25, arg1},
    {"C_multivariate_hmm", (DL_FUNC) &multivariate_hmm, 27, arg2},
    {"C_univariate_cleanup", (DL_FUNC) &univariate_cleanup, 0, NULL},
    {"C_multivariate_cleanup", (DL_FUNC) &multivariate_cleanup, 1, arg4},
    {"C_array3D_which_max", (DL_FUNC) &array3D_which_max, 3, arg5},
    {"C_array2D_mean", (DL_FUNC) &array2D_mean, 3, arg6},
    {"C_array3D_mean", (DL_FUNC) &array3D_mean, 3, arg7},
    {NULL, NULL, 0, NULL}
};


extern "C" {
void R_init_chromstaR(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
// 	R_forceSymbols(dll, TRUE);
}
}
