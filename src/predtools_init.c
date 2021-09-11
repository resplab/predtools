#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _predtools_Ccalc_mROC_stats(SEXP, SEXP);
extern SEXP _predtools_Csimulate_null_mROC_stats_unconditional(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_predtools_Ccalc_mROC_stats",                        (DL_FUNC) &_predtools_Ccalc_mROC_stats,                        2},
  {"_predtools_Csimulate_null_mROC_stats_unconditional", (DL_FUNC) &_predtools_Csimulate_null_mROC_stats_unconditional, 2},
  {NULL, NULL, 0}
};

void R_init_predtools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}