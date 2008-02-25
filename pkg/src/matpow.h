#include <R.h>
#include <Rinternals.h>

void matpow(double *x, int n, int k, double *z);

/* Use this as R - C API  -- so it's usable even from fortran user code : */
int F77_SUB(matpowf)(double *x, int *n, int *k, double *z);
