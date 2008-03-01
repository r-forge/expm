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
#include <R_ext/BLAS.h>

/* Workaround typo in R (<= 2.6.2)'s Lapack.h : *
 * in R >= 2.6.3, we just have
 #include <R_ext/Lapack.h>
*/
#include <R_ext/RS.h>
#define zlange_ __never__existing
#include <R_ext/Lapack.h>
#undef zlange_
extern double
F77_NAME(zlange)(const char *norm, const int *m, const int *n,
                 const Rcomplex *a, const int *lda,
                 double *work);
/* end{workaround} */

#include "locale.h"
#include "expm.h"

SEXP do_expm_eigen(SEXP x, SEXP tolin);
void expm_eigen(double *x, int n, double *z, double tol);

/* LAPACK functions *NOT* (yet) in R's Lapack library (nor Lapack.h),
 * see Lapack.h and html manual
 * <http://www.mathkeisan.com/UsersGuide/E/> : */

extern void
F77_NAME(zgecon)(const char* norm, const int* n,
                      const Rcomplex* a, const int* lda,
                      const double* anorm, double* rcond,
                      Rcomplex* work, double* rwork, int* info);
