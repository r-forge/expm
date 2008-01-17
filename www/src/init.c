/*
 *  Native routines registration
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "expm.h"

static const R_CallMethodDef CallEntries[] = {
    {"do_expm", (DL_FUNC) &do_expm, 1},
    {NULL, NULL, 0}
};

void R_init_expm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_RegisterCCallable("expm", "expm", (DL_FUNC) expm);
}
