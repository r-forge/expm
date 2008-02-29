/*  ===== File part of R package expm =====
 *
 *  expm-eigen.h  
 *
 *  Created by Christophe Dutang on 27/02/08.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "locale.h"
#include "expm.h"

SEXP do_expm_eigen(SEXP x, SEXP tolin);
void expm_eigen(double *x, int n, double *z, double tol);

/* Fortran auxiliary function, 
 * see BLAS.h or Lapack.h and html manual
 * <http://www.mathkeisan.com/UsersGuide/E/> */
extern void 
F77_NAME(zgemm)(const char *transa, const char *transb, const int *m,
                const int *n, const int *k, const Rcomplex *alpha,
                const Rcomplex *a, const int *lda,
                const Rcomplex *b, const int *ldb,
                const Rcomplex *beta, Rcomplex *c, const int *ldc);

extern void 
F77_NAME(zgecon)(const char* norm, const int* n,
                      const Rcomplex* a, const int* lda,
                      const double* anorm, double* rcond,
                      Rcomplex* work, double* rwork, int* info);
/*
extern double
F77_NAME(zlange)(const char *norm, const int *m, const int *n, 
                 const Rcomplex *a, const int *lda, 
                 double *work);*/