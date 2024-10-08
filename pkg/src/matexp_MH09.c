/* Copyright (C) 2014-2024 Martin Maechler
   Copyright (C) 2013-2014 Drew Schmidt

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, see <http://www.gnu.org/licenses/>.
*/


/* Matrix exponentiation algorithm from:
   "New Scaling and Squaring Algorithm for the Matrix Exponential", by
   Awad H. Al-Mohy and Nicholas J. Higham, August 2009
*/

#include <stdlib.h>
// #include <assert.h>

#include "expm.h"

#define SGNEXP(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))

// --------------------------------------------------------
// Utilities
// --------------------------------------------------------

// C = A * B for square matrices
static void matprod(int n, double *A, double *B, double *C)
{
    const double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, A, &n, B, &n, &zero, C, &n FCONE FCONE);
}



// Copy A ONTO B, i.e. B = A
static inline void matcopy(int n, double *A, double *B)
{
  F77_CALL(dlacpy)("A", &n, &n, A, &n, B, &n FCONE);
  //         uplo = ^
}



/** Identity matrix
 *
 * @param n  integer >= 1
 * @param a  n x n pre-allocated to contain the identity matrix
 */
static inline void mateye(const unsigned int n, double *a)
{
  int i;

  for (i=0; i<n*n; i++)
    a[i] = 0.0;

  i = 0;
  while (i < n*n)
  {
    a[i] = 1.0;
    i += n+1;
  }
}



// 1-norm for a square matrix
static double matnorm_1(const double *x, const int n)
{
  double norm = 0; // norm := max(colSums(abs(x)))
  for (int j=0; j<n; j++) {
      double tmp = 0;
      for (int i=0; i<n; i++)
	  tmp += fabs(x[i + j*n]);
      if (tmp > norm)
	  norm = tmp;
  }
  return norm;
}


#define NTHETA 5

static int matexp_scale_factor(const double *x, const int n)
{
    const double theta[] = {1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0};
    const double x_1 = matnorm_1(x, n);

    // MM: the following seems logically wrong: the theta[0:3] are *not* really relevant, just need
    for (int i=0; i < NTHETA; i++) {
	if (x_1 <= theta[i])
	    return 0;
    }

    int i = (int) ceil(log2(x_1/theta[4]));
    return 1 << i;
}

// ___ MM: FIXME  we have a  matpow() already in  ./matpow.c
//     --- Merge the two, keep the better one

// Matrix power by squaring: P = A^b (A is garbage on exit)
static void matpow_by_squaring(double *A, int n, int b, double *P)
{
    if (b == 1) {
	matcopy(n, A, P);
	return;
    }
    mateye(n, P);  // P := I
    if (b == 0)
	return;

    // General case: b >= 2
    double *TMP = (double *) R_alloc(n*n, sizeof(double));

    while (b) {
	if (b&1) { // P := P A
	    matprod(n, P, A, TMP);
	    matcopy(n, TMP, P);
	}

	b >>= 1;
	// A := A^2 :
	matprod(n, A, A, TMP);
	matcopy(n, TMP, A);
    }
}


// --------------------------------------------------------
// Matrix Exponentiation via Pade' Approximations
// --------------------------------------------------------

const double matexp_pade_coefs[14] =
{
  1.0,
  0.5,
  0.12,
  1.833333333333333333333e-2,
  1.992753623188405797101e-3,
  1.630434782608695652174e-4,
  1.035196687370600414079e-5,
  5.175983436853002070393e-7,
  2.043151356652500817261e-8,
  6.306022705717595115002e-10,
  1.483770048404140027059e-11,
  2.529153491597965955215e-13,
  2.810170546219962172461e-15,
  1.544049750670308885967e-17
};



/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j

   and q_m(x) = p_m(-x)
*/

// Workhorse for matexp_pade -- given C, do B := C & update (N, D)
void matexp_pade_fillmats(const int m, const int n, const int i,
			  double *N, double *D, double *B, const double *C)
{
  const double tmp = matexp_pade_coefs[i];
  const int sgn = SGNEXP(-1, i);

    /* Performs the following actions:
        B = C
        N +=          pade_coef[i] * C
        D += (-1)^i * pade_coef[i] * C
    */
    for (int j=0; j < m*n; j++) {
	double t_j = C[j]; B[j] = t_j;
	t_j *= tmp;
	N[j] +=     t_j;
	D[j] += sgn*t_j;
    }
}



/**
 * Exponentiation via Pade' expansion
 *
 * @param n
 * @param p
 * @param A
 * @param N
 */
static void matexp_pade(int n, const int p, double *A, double *N)
{
    int i, info = 0, n2 = n*n;
    // FIXME: check n2 (or n, such that n2 did not overflow !)

    // Power of A
    double *B = (double*) R_alloc(n2, sizeof(double));
    // Temporary storage for matrix multiplication;  matcopy(n, A, C);
    double *C = Memcpy((double*)R_alloc(n2, sizeof(double)), A, n2);
    double *D = (double*) R_alloc(n2, sizeof(double));

    for (i=0; i<n*n; i++) {
	N[i] = 0.0;
	D[i] = 0.0;
    }

    i = 0;
    while (i < n*n) {
	N[i] = 1.0;
	D[i] = 1.0;

	i += n+1;
    }


    // Fill N and D
    for (i=1; i<=p; i++)
    {
	// C = A*B
	if (i > 1)
	    matprod(n, A, B, C);

	// Update matrices
	matexp_pade_fillmats(n, n, i, N, D, B, C);
    }

    // return  N := inverse(D) %*% N -------------------
    int *ipiv = (int *) R_alloc(n, sizeof(int));
    /* assert(ipiv != NULL); */
    F77_CALL(dgesv)(&n, &n, D, &n, ipiv, N, &n, &info);

} // matexp_pade()


/**
 * Matrix Exponential
 *
 * @param x Input (square) matrix.  On exit, the values in x are "garbage"!
 * @param n Number of rows/cols of (square) matrix x.
 * @param p Order of the Pade' approximation. 0 < p <= 13.
 * @param ret On exit, ret = expm(x).
 */
void matexp_MH09(double *x, int n, const int p, double *ret)
{
  int m = matexp_scale_factor(x, n);

  if (m == 0) {
      matexp_pade(n, p, x, ret);
      return;
  }
  // else,  m >= 1 :

  int nn = n*n, one = 1;
  double tmp = 1. / ((double) m);

  F77_CALL(dscal)(&nn, &tmp, x, &one);

  matexp_pade(n, p, x, ret);

  matcopy(n, ret, x);

  matpow_by_squaring(x, n, m, ret);
}

#if 0 // >>>>>>>>>>> FIXME: Not yet, does not even compile  -- (and hence not yet called) <<<<<<<<<

void matexp_MH09_z(Rcomplex *x, int n, const int p, Rcomplex *ret)
{
  int m = matexp_scale_factor(x, n);

  if (m == 0) {
      matexp_pade(n, p, x, ret);
      return;
  }

  int nn = n*n, one = 1;
  Rcomplex tmp = 1. / ((Rcomplex) m);

  F77_CALL(dscal)(&nn, &tmp, x, &one);

  matexp_pade(n, p, x, ret);
  matcopy(n, ret, x);
  matpow_by_squaring(x, n, m, ret);
}
#endif


// --------------------------------------------------------
// R Wrapper
// --------------------------------------------------------

SEXP R_matexp_MH09(SEXP x, SEXP p)
{
    const int n = nrows(x), n2 = n*n;
    SEXP x_ = duplicate(x); PROTECT_INDEX xpi;
    PROTECT_WITH_INDEX(x_, &xpi);
    Rboolean is_C = isComplex(x);
    if(is_C) {
    } else if(!isReal(x)) /* coercion to 'Real' numeric */
	REPROTECT(x_ = coerceVector(x_, REALSXP), xpi);
    SEXP R = PROTECT(allocMatrix(TYPEOF(x_), n, n));
    if(isComplex(x)) {
	error(_("matexp_MH09(.) is _not yet_ implemented for complex matrices"));
#if 0 // see above
	// MM{FIXME}: We have already duplicated x to x_ : this should *not* be needed:
	Rcomplex *x_cp = Memcpy((Rcomplex*)R_alloc(n2, sizeof(Rcomplex)), COMPLEX(x_), n2);
	matexp_MH09_z(x_cp, n, INTEGER(p)[0], COMPLEX(R));
# endif

    } else { // isReal(x_)
	// MM{FIXME}: We have already duplicated x to x_ : this should *not* be needed:
	double *x_cp = Memcpy((double*)R_alloc(n2, sizeof(double)), REAL(x_), n2);
	matexp_MH09(x_cp, n, INTEGER(p)[0], REAL(R));
    }
    UNPROTECT(2);
    return R;
}

