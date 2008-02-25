
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "locale.h"

SEXP do_expm(SEXP x);
void expm(double *x, int n, double *z);

