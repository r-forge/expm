/*
 *  Native routines registration
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "expm.h"
#include "matpow.h"

static const R_CallMethodDef CallEntries[] = {
    {"do_expm", (DL_FUNC) &do_expm, 1},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
    {"matpowf", (DL_FUNC) &F77_SUB(matpowf), 4},
    {NULL, NULL, 0}
};

void R_init_expm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortEntries, NULL);
    /* callable C code from other packages C code :*/
    R_RegisterCCallable("expm", "expm", (DL_FUNC) expm);
    R_RegisterCCallable("matpow", "matpow", (DL_FUNC) matpow);
}
