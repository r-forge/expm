/*  ===== File part of R package expm =====
 *
 *  Function to compute the matrix exponential
 *
 *     exp(M) = sum(n = 0:Inf; M^n / n!),
 *
 *  where M is an (n x n) matrix.
 *
 *  The functions therein use LAPACK and BLAS routines. Nicely
 *  formatted man pages for these can be found at
 *
 *    <http://www.mathkeisan.com/UsersGuide/E/>
 *
 *  AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Christophe
 *  Dutang, based on code in package Matrix.
 */

#include "expm.h"


/* For matrix exponential calculations. Pade constants
 *
 *   n_{pqj} = [(p + q - j)! p!]/[(p + q)! j! (p - j)!]
 *
 * and
 *
 *   d_{pqj} = [(p + q - j)! q!]/[(p + q)! j! (q - j)!]
 *
 * for p = q = 8 and j = 1, ..., 8.
 */
const static double padec88 [] =
{
    5.0000000000000000e-1,
    1.1666666666666667e-1,
    1.6666666666666667e-2,
    1.6025641025641026e-3,
    1.0683760683760684e-4,
    4.8562548562548563e-6,
    1.3875013875013875e-7,
    1.9270852604185938e-9,
};


/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. Mostly lifted from the core of fonction
 * expm() of package Matrix, which is itself based on the function of
 * the same name in Octave.
 */
void expm(double *x, int n, double *z)
{
    if (n == 1)
	z[0] = exp(x[0]);		/* scalar exponential */
    else
    {
	/* Constants */
	const double one = 1.0, zero = 0.0;
	const int i1 = 1, nsqr = n * n, np1 = n + 1;
	/* Variables */
	int i, j, is_uppertri = TRUE;;
	int ilo, ihi, iloscal, ihiscal, info, sqrpowscal;
	double infnorm, trshift, m1pj = -1;

	/* Arrays */
	int *pivot    = (int *) R_alloc(n, sizeof(int)); /* pivot vector */
	double *perm  = (double *) R_alloc(n, sizeof(double)); /* permutation array */
	double *scale = (double *) R_alloc(n, sizeof(double)); /* scale array */
	double *work  = (double *) R_alloc(nsqr, sizeof(double)); /* workspace array */
	double *npp   = (double *) R_alloc(nsqr, sizeof(double)); /* num. power Pade */
	double *dpp   = (double *) R_alloc(nsqr, sizeof(double)); /* denom. power Pade */

	Memcpy(z, x, nsqr);

	/* Check if matrix x is upper triangular; stop checking as
	 * soon as a non-zero value is found below the diagonal. */
	for (i = 0; i < n - 1 && is_uppertri; i++)
	    for (j = i + 1; j < n; j++)
		if (!(is_uppertri = x[i * n + j] == 0.0))
		    break;


	/* Step 1 of preconditioning: shift diagonal by average diagonal. */
	trshift = 0.0;
	for (i = 0; i < n; i++)
	    trshift += x[i * np1];
	trshift /= n;		/* average diagonal element */
	if (trshift > 0.0)
	    for (i = 0; i < n; i++)
		z[i * np1] -= trshift;

	/* Step 2 of preconditioning: balancing with dgebal. */
	if (is_uppertri) {
	    /* no need to permute if x is upper triangular */
	    ilo = 1;
	    ihi = n;
	}
	else {
	    F77_CALL(dgebal)("P", &n, z, &n, &ilo, &ihi, perm, &info);
	    if (info)
		error(_("LAPACK routine dgebal returned info code %d when permuting"), info);
	}
	F77_CALL(dgebal)("S", &n, z, &n, &iloscal, &ihiscal, scale, &info);
	if (info)
	    error(_("LAPACK routine dgebal returned info code %d when scaling"), info);

	/* Step 3 of preconditioning: Scaling according to infinity
	 * norm (a priori always needed). */
	infnorm = F77_CALL(dlange)("I", &n, &n, z, &n, work);
	sqrpowscal = (infnorm > 0) ? imax2((int) 1 + log(infnorm)/log(2.0), 0) : 0;
	if (sqrpowscal > 0) {
	    double scalefactor = R_pow_di(2, sqrpowscal);
	    for (i = 0; i < nsqr; i++)
		z[i] /= scalefactor;
	}

	/* Pade approximation (p = q = 8): compute x^8, x^7, x^6,
	 * ..., x^1 */
	for (i = 0; i < nsqr; i++)
	{
	    npp[i] = 0.0;
	    dpp[i] = 0.0;
	}
	for (j = 7; j >= 0; j--)
	{
	    /* npp = z * npp + padec88[j] * z */
	    F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, npp,
			     &n, &zero, work, &n);
	    /* npp <- work + padec88[j] * z */
	    for (i = 0; i < nsqr; i++)
		npp[i] = work[i] + padec88[j] * z[i];
	    /* dpp = z * dpp + (-1)^j * padec88[j] * z */
	    F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, dpp,
			     &n, &zero, work, &n);
	    for (i = 0; i < nsqr; i++)
		dpp[i] = work[i] + m1pj * padec88[j] * z[i];
	    m1pj *= -1; /* (-1)^j */
	}
	/* power 0 */
	for (i = 0; i < nsqr; i++)
	    dpp[i] *= -1.0;
	for (j = 0; j < n; j++)
	{
	    npp[j * np1] += 1.0;
	    dpp[j * np1] += 1.0;
	}

	/* Pade approximation is (dpp)^-1 * npp. */
	F77_CALL(dgetrf) (&n, &n, dpp, &n, pivot, &info);
	if (info)
	    error(_("LAPACK routine dgetrf returned info code %d"), info);
	F77_CALL(dgetrs) ("N", &n, &n, dpp, &n, pivot, npp, &n, &info);
	if (info)
	    error(_("LAPACK routine dgetrs returned info code %d"), info);

	Memcpy(z, npp, nsqr);

	/* Now undo all of the preconditioning */
	/* Preconditioning 3: square the result for every power of 2 */
	while (sqrpowscal--)
	{
	    F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, z, &n,
			    z, &n, &zero, work, &n);
	    Memcpy(z, work, nsqr);
	}


	/* Preconditioning 2: Inversion of 'dgebal()' :
	 * ------------------ Note that dgebak() seems *not* applicable */

	/* Step 2 a)  apply inverse scaling */
	for (j = 0; j < n; j++)
	    for (i = 0; i < n; i++)
		z[i + j * n] *= scale[i]/scale[j];

	/* 2 b) Inverse permutation  (if not the identity permutation) */

	if (ilo != 1 || ihi != n) {

	    /* ---- new code by Martin Maechler ----- */

#define SWAP_ROW(I,J) F77_CALL(dswap)(&n, &z[(I)], &n, &z[(J)], &n)

#define SWAP_COL(I,J) F77_CALL(dswap)(&n, &z[(I)*n], &i1, &z[(J)*n], &i1)

#define RE_PERMUTE(I)				\
	    int p_I = (int) (perm[I]) - 1;	\
	    SWAP_COL(I, p_I);			\
	    SWAP_ROW(I, p_I)

	    /* reversion of "leading permutations" : in reverse order */
	    for (i = (ilo - 1) - 1; i >= 0; i--) {
		RE_PERMUTE(i);
	    }

	    /* reversion of "trailing permutations" : applied in forward order */
	    for (i = (ihi + 1) - 1; i < n; i++) {
		RE_PERMUTE(i);
	    }
	}

	/* Preconditioning 1: Trace normalization */
	if (trshift > 0)
	{
	    double mult = exp(trshift);
	    for (i = 0; i < nsqr; i++)
		z[i] *= mult;
	}

    }
}

/* Main function, the only one used by .External(). */
SEXP do_expm(SEXP x)
{
    SEXP dims, z;
    int n;
    double *rx = REAL(x), *rz;

    if (!isNumeric(x) || !isMatrix(x))
	error(_("invalid argument: not a numeric matrix"));

    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    if (n != INTEGER(dims)[1])
	error(_("non-square matrix"));
    if (n == 0)
	return(allocMatrix(REALSXP, 0, 0));

    PROTECT(z = allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    expm(rx, n, rz);

    setAttrib(z, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);
    return z;
}
