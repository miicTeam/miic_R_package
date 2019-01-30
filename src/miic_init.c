#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mydiscretizeMDL(SEXP, SEXP);
extern SEXP mydiscretizeMutual(SEXP, SEXP, SEXP, SEXP, SEXP,
                               SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP evaluateEffn(SEXP, SEXP, SEXP);
extern SEXP orientationProbability(SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP skeleton(SEXP, SEXP, SEXP, SEXP, SEXP, 
                     SEXP, SEXP, SEXP, SEXP, SEXP, 
                     SEXP, SEXP, SEXP, SEXP, SEXP, 
                     SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mydiscretizeMDL",        (DL_FUNC) &mydiscretizeMDL,          2},
    {"mydiscretizeMutual",     (DL_FUNC) &mydiscretizeMutual,       10},
    {"evaluateEffn",           (DL_FUNC) &evaluateEffn,             3},
    {"orientationProbability", (DL_FUNC) &orientationProbability,  15},
    {"skeleton",               (DL_FUNC) &skeleton,                20},
    {NULL, NULL, 0}
};

void R_init_miic(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
