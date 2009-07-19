#include <R.h>
#include <Rdefines.h>
#include "newuoa_head.h"

OptStruct OS;

SEXP newuoa_c(SEXP par_arg, SEXP fn, SEXP control, SEXP rho)
{
    int     i, j;
    int    n, npt, maxfun, iprint, wsize;
    double rhobeg, rhoend;

    double  *par;
    
    SEXP    sexp_rss, sexp_feval;
    SEXP    out, out_names;
    SEXP env;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));

    OS->feval = 0; 

    PROTECT(OS->par = duplicate(par_arg));
    n = length(OS->par);

    switch (TYPEOF(OS->par)) {
    case REALSXP:
      break;
    default:
      error("`par' that you provided is non-list and non-numeric!");
    }
    
    if (!isFunction(fn)) error("fn is not a function!");
    PROTECT(OS->fcall = lang2(fn, OS->par));

    if (!isEnvironment(rho)) error("rho is not an environment!");
    OS->env = rho;

    par         = real_vector(n);
    if (IS_NUMERIC(OS->par)) {
      for (i = 0; i < n; i++)
	par[i] = NUMERIC_POINTER(OS->par)[i];
    }
    rhobeg    = NUMERIC_VALUE(getListElement(control, "rhobeg"));
    rhoend    = NUMERIC_VALUE(getListElement(control, "rhoend"));

    npt    = NUMERIC_VALUE(getListElement(control, "npt"));
    maxfun    = NUMERIC_VALUE(getListElement(control, "maxfun"));
    iprint    = NUMERIC_VALUE(getListElement(control, "iprint"));
    wsize    = NUMERIC_VALUE(getListElement(control, "wsize"));
    Rprintf("w \n",wsize);
    
    /*size w needs to be >= (((OS->npt+13) * (OS->npt+n) + (3*n*(n+3)/2)))*/
    double w[wsize]; 

/*========================================================================*/
    
    F77_CALL(newuoa)(&n, &npt, par, &rhobeg, &rhoend, 
		     &iprint, &maxfun, w);
    
/*========================================================================*/
    
    
    PROTECT(sexp_rss = NEW_NUMERIC(1));
    NUMERIC_POINTER(sexp_rss)[0] = OS->rss;
    
    PROTECT(sexp_feval = NEW_NUMERIC(1));
    NUMERIC_POINTER(sexp_feval)[0] = OS->feval;

    PROTECT(out = NEW_LIST(3));
    SET_VECTOR_ELT(out, 0, OS->par);
    SET_VECTOR_ELT(out, 1, sexp_rss);
    SET_VECTOR_ELT(out, 2, sexp_feval);

    PROTECT(out_names = NEW_STRING(3));
    SET_STRING_ELT(out_names, 0, mkChar("par"));
    SET_STRING_ELT(out_names, 1, mkChar("fval"));
    SET_STRING_ELT(out_names, 2, mkChar("feval"));
    
    SET_NAMES(out, out_names);

    UNPROTECT(6);

    return out;
}
 
