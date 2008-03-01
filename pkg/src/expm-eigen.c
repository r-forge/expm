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
 *  AUTHORS: Christophe Dutang, based on code eigen,
 *
 *  i.e., function 'modLa_rg' and 'modLa_dgesv' in R's
 *  <Rsource>/src/modules/lapack/lapack.c, used in eigen()
 */



#include "expm-eigen.h"

void expm_eigen(double *x, int n, double *z, double tol)
{
    if (n == 1)
        z[0] = exp(x[0]);		/* scalar exponential */
    else
    {
        /* Constants */
        int i, j;
        int nsqr = n * n, np1 = n + 1;
        int info, lwork, is_conjug, is_diag;
        double onenorm, rcond, tmp;
        char jobVL[1], jobVR[1];
        char *transa = "N";
        char *one = "1";
        Rcomplex cone, czero;

        /* Arrays */
        int *ipiv = (int *) R_alloc(n, sizeof(int)); /* permutation vector */
        double *left, *right, *workdiag; /* left and right eigenvectors and workspace for diagonalisation */
        double *wR = (double *) R_alloc(n, sizeof(double)); /* real part of eigenvalues */
        double *wI = (double *) R_alloc(n, sizeof(double)); /* imaginary part of eigenvalues */
        double *rworksing = (double *) R_alloc(2*n, sizeof(double)); /* working vector to test the singularity */
        Rcomplex *eigvect = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* (right) eigenvectors matrix */
        Rcomplex *eigvectinv = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* its inverse */
        Rcomplex *expeigval; /* complex matrix diag(exp(eigenvalues)) */
        Rcomplex *ctmp = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* temp working variable */
        Rcomplex *worksing = (Rcomplex *) R_alloc(2*n, sizeof(Rcomplex)); /* workspace to test the singularity */

        R_CheckStack();

        Memcpy(z, x, nsqr);

        /* Test if x is diagonalisable by computing its eigenvalues and (right) eigenvectors */
        /* code based on modLa_rg in lapack.c, used in eigen.R */
        jobVL[0] = 'N';
        left = (double *) 0;
        jobVR[0] = 'V';
        right = (double *) R_alloc(nsqr, sizeof(double));

        /* 1 - ask for optimal size of work array */
        lwork = -1;
        F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI, left, &n, right, &n, &tmp, &lwork, &info);
        if (info != 0)
            error(_("error code %d from Lapack routine dgeev"), info);
        lwork = (int) tmp;
        workdiag = (double *) R_alloc(lwork, sizeof(double));

        /* 2 - compute eigenvalues and (right) eigenvectors */
        F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI, left, &n, right, &n, workdiag, &lwork, &info);
        if (info != 0)
            error(_("error code %d from Lapack routine dgeev"), info);

        /*Rprintf("eigen value\n");
        for(i = 0; i < n; i++)
            Rprintf("%f +i* %f\n", wR[i], wI[i]);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                Rprintf("%f\n", right[i * n + j]);*/

        /* try to invert the eigenvectors matrix */
        /* 1 - build the Rcomplex matrix with eigenvectors */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                is_conjug = 0;
                if(i < n-1)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i+1] && wI[i] == -wI[i+1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[i * n + j];
                        eigvect[i * n + j].i = right[(i+1) * n + j];
                    }
                }
                if(i > 0)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i-1] && wI[i] == -wI[i-1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[(i-1) * n + j];
                        eigvect[i * n + j].i = -right[i * n + j];
                    }
                }
                /* real eigenvalues */
                        //Rprintf("%d-%d-%d\n",i,j,is_conjug);
                if(!is_conjug)
                {
                    eigvect[i * n + j].r = right[i * n + j];
                    eigvect[i * n + j].i = 0.0;
                }
                /* eigvectinv initialise with the identity matrix */
                eigvectinv[i * n +j].r = (i == j) ? 1.0 : 0.0;
                eigvectinv[i * n +j].i = 0.0;
            }
        }

        /*Rprintf("eigvect\n");
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                Rprintf("%.10f +i*%.10f\n", eigvect[i * n + j].r, eigvect[i * n + j].i);
        Rprintf("eigvectinv\n");
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                Rprintf("%.10f +i*%.10f\n", eigvectinv[i * n + j].r, eigvectinv[i * n + j].i);*/

        /* 2 - store the matrix eigvect (because function zgesv will change it) */
        Memcpy(ctmp, eigvect, nsqr);

        /* 3 - solve a linear complex equation system with eigvectinv equals
         * to matrix identity. hence, on exit eigvectinv contains the
         * inverse of complex matrix eigvect. code base on solve.R */
        F77_CALL(zgesv)(&n, &n, eigvect, &n, ipiv, eigvectinv, &n, &info);

        if (info > 0)
            is_diag = 0; //matrix eigvect is exactly singular.
        if (info < 0)
            error(_("argument %d of Lapack routine dgesv had invalid value"), -info);
        if (info == 0)
            is_diag = 1;

        /* check if matrix eigvect is numerically singular */
        /* 1 - compute the one norm of the complex matrix eigvect */
        onenorm = F77_CALL(zlange)("1", &n, &n, eigvect, &n, (double*) NULL);

        /* 2 - estimates the reciprocal of the condition number of the one norm of eigvect */

	/* This is *NOT* part of R's LAPACK : */

/*         F77_CALL(zgecon)("1", &n, eigvect, &n, &onenorm, &rcond, */
/* 			 worksing, rworksing, &info); */

/*         if (rcond < tol) */
/*         { */
/*             Rprintf("computationally singular\n"); */
/*             is_diag=0; */
/*         } */

        if (is_diag)
        {
            /* x is diagonalisable so
             * compute complex matrix operations :
             * eigvect %*% diag(exp(eigenvalues)) %*% eigvectinv */

            /* 1 - expeigval is the complex matrix diag(exp(eigenvalues)) */
            /* code based on z_exp in complex.c */
            expeigval = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex));
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if(i == j)
                    {
                        expeigval[i * n +j].r = exp(wR[i]) * cos(wI[i]);
                        expeigval[i * n +j].i = exp(wR[i]) * sin(wI[i]);
                    }
                    else
                    {
                        expeigval[i * n +j].r = 0.0;
                        expeigval[i * n +j].i = 0.0;
                    }
                }
            }

            /*Rprintf("expeigval\n");
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Rprintf("%.10f +i*%.10f\n", expeigval[i * n + j].r, expeigval[i * n + j].i);*/

            /* 2 - restore the matrix eigvect */
            Memcpy(eigvect, ctmp, nsqr);

            /*
            Rprintf("\nAPRES\n\n");

            Rprintf("eigvect\n");
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Rprintf("%.10f +i*%.10f\n", eigvect[i * n + j].r, eigvect[i * n + j].i);
            Rprintf("eigvectinv\n");
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Rprintf("%.10f +i*%.10f\n", eigvectinv[i * n + j].r, eigvectinv[i * n + j].i);*/

            cone.r = 1.0; cone.i = czero.r = czero.i = 0.0;
            /* 3 - compute (complex) matrix product: ctmp <- eigvect * expeigval */
            F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, eigvect, &n, expeigval, &n, &czero, ctmp, &n);
            /*
            Rprintf("eigvect * expeigval\n");
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Rprintf("%.10f +i*%.10f\n", ctmp[i * n + j].r, ctmp[i * n + j].i);*/

            /* 4 - compute (complex) matrix product: expeigval <- ctmp * eigvectinv */
            F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, ctmp, &n, eigvectinv, &n, &czero, expeigval, &n);
            /*
            Rprintf(" ctmp * eigvectinv\n");
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Rprintf("%.10f +i*%.10f\n", expeigval[i * n + j].r, expeigval[i * n + j].i);*/


            /* store the real part in z */
            /* the matrix exponential is real, even if x has complex eigen values. */
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    z[i * n + j] = expeigval[i * n + j].r;

        }
	else
	    expm(x, n, z, Ward_2);
    }
}

/* Main function, the only one used by .Call(). */
SEXP do_expm_eigen(SEXP x, SEXP tolin)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;
    double tol = asReal(tolin);

    if (!isNumeric(x) || !isMatrix(x))
        error(_("invalid argument"));

    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m)
        error(_("non-square matrix"));
    if (n == 0)
        return(allocVector(REALSXP, 0));

    PROTECT(z = allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    expm_eigen(rx, n, rz, tol);

    setAttrib(z, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);

    return z;
}
