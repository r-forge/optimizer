#include <R.h>
#include <Rdefines.h>
#include "bobyqa_head.h"

OptStruct OS;

SEXP bobyqa_c(SEXP par_arg, SEXP xl_arg, SEXP xu_arg, SEXP fn, SEXP control, 
	      SEXP rho)
{
    int     i, j;
    int    n, npt, maxfun, iprint, wsize;
    double rhobeg, rhoend;

    double  *par, *xl, *xu;
    
    SEXP    sexp_rss, sexp_feval;
    SEXP    out, out_names;
    SEXP env;


    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));

    OS->feval = 0; 

    PROTECT(OS->par = duplicate(par_arg));

    n = length(OS->par);

    if(n != length(xl_arg) || n != length(xu_arg))
      error("parameter bounds provided have the wrong length!");
    
    switch (TYPEOF(OS->par)) {
    case REALSXP:
      break;
    default:
      error("`par' that you provided is non-list and non-numeric!");
    }
    switch (TYPEOF(xl_arg)) {
    case REALSXP:
      break;
    default:
      error("`xl' that you provided is non-list and non-numeric!");
    }
    switch (TYPEOF(xu_arg)) {
    case REALSXP:
      break;
    default:
      error("`xu' that you provided is non-list and non-numeric!");
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
    xl         = real_vector(n);
    if (IS_NUMERIC(xl_arg)) {
      for (i = 0; i < n; i++)
	xl[i] = NUMERIC_POINTER(xl_arg)[i];
    }
    xu         = real_vector(n);
    if (IS_NUMERIC(xu_arg)) {
      for (i = 0; i < n; i++)
	xu[i] = NUMERIC_POINTER(xu_arg)[i];
    }

    rhobeg    = NUMERIC_VALUE(getListElement(control, "rhobeg"));
    rhoend    = NUMERIC_VALUE(getListElement(control, "rhoend"));

    npt    = NUMERIC_VALUE(getListElement(control, "npt"));
    maxfun    = NUMERIC_VALUE(getListElement(control, "maxfun"));
    iprint    = NUMERIC_VALUE(getListElement(control, "iprint"));
    wsize    = NUMERIC_VALUE(getListElement(control, "wsize"));
    
    /*size w needs to be >= (NPT+5)*(NPT+N)+3*N*(N+5)/2. */
    double w[wsize]; 

/*========================================================================*/
    
    F77_CALL(bobyqa)(&n, &npt, par, xl, xu, &rhobeg, &rhoend, 
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
 
