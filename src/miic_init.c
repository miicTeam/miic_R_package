#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h>  // for NULL

// FIXME: Check these declarations against the C/Fortran source code.
// .Call calls
extern SEXP mydiscretizeMDL(SEXP, SEXP);
extern SEXP mydiscretizeMutual(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP evaluateEffn(SEXP, SEXP, SEXP);
extern SEXP reconstruct(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mydiscretizeMDL", (DL_FUNC)&mydiscretizeMDL, 2},
    {"mydiscretizeMutual", (DL_FUNC)&mydiscretizeMutual, 11},
    {"evaluateEffn", (DL_FUNC)&evaluateEffn, 3},
    {"reconstruct", (DL_FUNC)&reconstruct, 25}, {NULL, NULL, 0}};

void R_init_miic(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
