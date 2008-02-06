#include <R.h>
#include <Rinternals.h>

SEXP do_expm(SEXP x);
void expm(double *x, int n, double *z);

extern void
F77_NAME(zgemm)(const char *transa, const char *transb, const int *m,
	    const int *n, const int *k, const Rcomplex *alpha,
	    const Rcomplex *a, const int *lda,
	    const Rcomplex *b, const int *ldb,
	    const Rcomplex *beta, Rcomplex *c, const int *ldc);