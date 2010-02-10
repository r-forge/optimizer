#include <R.h>
#include <Rdefines.h>
#include "newuoa_head.h"

void F77_SUB(resfun)(int *n, double *par, double *f) {
  int i;
  SEXP ff;
  for (i = 0; i < *n; i++) {
    if (!R_FINITE(par[i]))
      error("non-finite parameter to user function");
    NUMERIC_POINTER(OS->par)[i] = par[i];
  }
  SETCADR(OS->fcall, OS->par);
  PROTECT(ff = eval(OS->fcall, OS->env));
  f[0] = NUMERIC_POINTER(ff)[0];
  OS->rss = f[0];
  if (!R_FINITE(OS->rss))
      error("non-finite function value from user function");
  OS->feval = OS->feval++;

  UNPROTECT(1);	
  
}
