/* --- this is basically a copy from  the code in  CRAN package actuar's
 * src/util.c  (by Vincent Goulet)
 * with slight shortcuts by Martin Maechler: */

/* Power of a matrix x^k := x x ... x, where x in an (n x n) matrix
 * and k is an *integer* (including -1). This function is fairly naive
 * with little error checking since it is currently used in a very
 * narrow and already controlled context.
 */
#include <R_ext/BLAS.h>

#include "matpow.h"

void matpow(double *x, int n, int k, double *z);

/* Use this as R - C API  -- so it's usable even from fortran user code : */

/* TODO : replace this with .Call() version !
 * void R_matpow(double *x, int *n, int *k, double *z) */

int F77_SUB(matpowf)(double *x, int *n, int *k, double *z)
{
    int kk = *k; /* pass an 'own' integer; it will be modified below! */
    matpow(x, *n, kk, z);
    return 0;
}

void matpow(double *x, int n, int k, double *z)
{
    if (k == 0) { /* Return identity matrix */
	int i, j;
	for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
		z[i * n + j] = (i == j) ? 1.0 : 0.0;
	return;
    }
    else if (k < 0) {

	REprintf("If you *REALLY* need <matrix> ^ (-n), use (m^n)^(-1) yourself!");
    }
    else { /* k >= 1 */

	static const char *transa = "N";
	static const double one = 1.0, zero = 0.0;
	int nSqr = n * n;
	double /* temporary matrix */
	    *tmp  = (double *) R_alloc(nSqr, sizeof(double));

	/* Take powers in multiples of 2 until there is only one
	 * product left to make. That is, if k = 5, compute (x * x),
	 * then ((x * x) * (x * x)) and finally ((x * x) * (x * x)) *
	 * x. Idea taken from Octave in file .../src/xpow.cc. */
	Memcpy(z, x, (size_t) nSqr);

	k--;
	while (k > 0) {
	    if (k & 1) {	/* z := z * x */
		F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
				z, &n, x, &n, &zero, tmp, &n);
		Memcpy(z, tmp, (size_t) nSqr);
	    }

	    if(k == 1)
		break;
	    k >>= 1; /* efficient division by 2 ; now have k >= 1 */
	    /* x := x * x */
	    F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
			    x, &n, x, &n, &zero, tmp, &n);
	    Memcpy(x, tmp, (size_t) nSqr);
	}
    }
}
