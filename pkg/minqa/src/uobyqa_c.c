#include <R.h>
#include <Rdefines.h>
#include "uobyqa_head.h"

/* Added 1 by D Bates to move output from Powell Fortran codes */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("minqa", String)
#else
#define _(String) (String)
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif
/* End 1 D Bates */



OptStruct OS;

SEXP uobyqa_c(SEXP par_arg, SEXP fn, SEXP control, SEXP rho)
{
    int     i, j;
    int    n, maxfun, iprint, wsize;
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

    maxfun    = NUMERIC_VALUE(getListElement(control, "maxfun"));
    iprint    = NUMERIC_VALUE(getListElement(control, "iprint"));
    wsize    = NUMERIC_VALUE(getListElement(control, "wsize"));
    
    /* must be > ( N**4 + 8*N**3 + 23*N**2 + 42*N + 
       max [ 2*N**2 + 4, 18*N ] ) / 4. */
    double w[wsize]; 

/*========================================================================*/
    
    F77_CALL(uobyqa)(&n, par, &rhobeg, &rhoend, 
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
 
/* Added 2 by JN following D Bates to move output from Powell Fortran codes

   The routines minqer, minquit, minqi3 and minqir already defined in bobyqa_c.c

 */
/* End JN following 2 D Bates */
/* Added J Nash */
 
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("minqa", String)
#else
#define _(String) (String)
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

void attribute_hidden
F77_NAME(minqiru)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	int i, nn = *n;
	Rprintf(_("At return from uobyqa\n"));
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}
/* End JN */


