#include <R.h>
#include <Rdefines.h>
#include "bobyqa_head.h"

void F77_SUB(resfunbobyqa)(int *n, double *par, double *f) {
  int i;
  SEXP ff;
  for (i = 0; i < *n; i++) {
    if (!R_FINITE(par[i]))
      error("non-finite value supplied by lmdif!");
    NUMERIC_POINTER(OS->par)[i] = par[i];
  }
  SETCADR(OS->fcall, OS->par);
  PROTECT(ff = eval(OS->fcall, OS->env));
  f[0] = NUMERIC_POINTER(ff)[0];
  OS->rss = f[0];
  OS->feval = OS->feval++;
  Rprintf("here\n\n\n");
  UNPROTECT(1);	
  
}
