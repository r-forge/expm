/* --- this is basically a copy from  the code in  CRAN package actuar's
 * src/util.c  (by Vincent Goulet)
 * with slight shortcuts by Martin Maechler: */

/* Power of a matrix x^k := x x ... x, where x in an (n x n) matrix
 * and k is an *integer* (including -1). This function is fairly naive
 * with little error checking since it is currently used in a very
 * narrow and already controlled context.
 */

#include "matpow.h"

/* .Call() this from R : */
SEXP R_matpow(SEXP x, SEXP k)
{
    if(!isMatrix(x))
	error(_("not a matrix"));
    else {
	SEXP dims = getAttrib(x, R_DimSymbol), z;
	int n = INTEGER(dims)[0],
	    kk = INTEGER(k)[0]; /* need copy, as it is altered in matpow() */

	if(!isNumeric(x))
	    PROTECT(x = coerceVector(x, REALSXP)); /* may give error,...*/
	else
	    PROTECT(x = duplicate(x)); /* since matpow() will alter it */

	if (n != INTEGER(dims)[1]) {
	    UNPROTECT(1);
	    error(_("non-square matrix"));
	}
	if (n == 0)
	    return(allocMatrix(REALSXP, 0, 0));

	PROTECT(z = allocMatrix(REALSXP, n, n));
	setAttrib(z, R_DimNamesSymbol,
		  getAttrib(x, R_DimNamesSymbol));

	matpow(REAL(x), n, kk, REAL(z));

	UNPROTECT(2);
	return z;
    }
}

void matpow(double *x, int n, int k, double *z)
/* Compute  z :=  x %^% k ;  x an (n x n) square "matrix" in column-order;
 * NB: x[] will be altered !! --- the caller must make a copy if needed */
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
