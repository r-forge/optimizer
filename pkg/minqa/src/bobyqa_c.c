#include <R.h>
#include <Rdefines.h>
#include "bobyqa_head.h"

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

SEXP bobyqa_c(SEXP par_arg, SEXP xl_arg, SEXP xu_arg, SEXP fn, SEXP control, 
	      SEXP rho)
{
    int i, j;
    int n = LENGTH(par_arg), npt, maxfun, iprint, wsize;
    double rhobeg, rhoend;
    double  *par, *xl, *xu;
    
    SEXP sexp_rss, sexp_feval, sexp_fval, out, out_names, env;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
    OS->feval = 0; 
    PROTECT(OS->par = coerceVector(par_arg, REALSXP));

    if(n != LENGTH(xl_arg) || n != LENGTH(xu_arg) ||
       !isReal(xl_arg) || !isReal(xu_arg))
	error(_("parameter bounds provided have the wrong length or mode"));
    
    if (!isFunction(fn)) error(_("fn is not a function!"));
    PROTECT(OS->fcall = lang2(fn, OS->par));

    if (!isEnvironment(rho)) error("rho is not an environment!");
    OS->env = rho;

    /* it is easier to extract by name from an environment than a list */
    if (!isEnvironment(control))
	error(_("control argument should be an environment"));
    rhobeg = asReal(findVarInFrame(control, install("rhobeg")));
    rhoend = asReal(findVarInFrame(control, install("rhoend")));
    npt    = asInteger(findVarInFrame(control, install("npt")));
    maxfun = asInteger(findVarInFrame(control, install("maxfun")));
    iprint = asInteger(findVarInFrame(control, install("iprint")));
    wsize  = asInteger(findVarInFrame(control, install("wsize")));
    
    /*size w needs to be >= (NPT+5)*(NPT+N)+3*N*(N+5)/2. */
    double *w = (double*) R_alloc(wsize, sizeof(double));
    double *fval = (double*) R_alloc(npt, sizeof(double)); 
    
/*========================================================================*/
    
    F77_CALL(bobyqa)(&n, &npt, REAL(OS->par), REAL(xl_arg), REAL(xu_arg),
		     &rhobeg, &rhoend, &iprint, &maxfun, w, fval);
    
/*========================================================================*/

 
    PROTECT(sexp_fval = NEW_NUMERIC(npt));
    for (i = 0; i < npt; i++)
      NUMERIC_POINTER(sexp_fval)[i] = fval[i];
 
    PROTECT(sexp_rss = NEW_NUMERIC(1));
    NUMERIC_POINTER(sexp_rss)[0] = OS->rss;
    
    PROTECT(sexp_feval = NEW_NUMERIC(1));
    NUMERIC_POINTER(sexp_feval)[0] = OS->feval;

    PROTECT(out = NEW_LIST(4));
    SET_VECTOR_ELT(out, 0, OS->par);
    SET_VECTOR_ELT(out, 1, sexp_rss);
    SET_VECTOR_ELT(out, 2, sexp_feval);
    SET_VECTOR_ELT(out, 3, sexp_fval);

    PROTECT(out_names = NEW_STRING(4));
    SET_STRING_ELT(out_names, 0, mkChar("par"));
    SET_STRING_ELT(out_names, 1, mkChar("fval"));
    SET_STRING_ELT(out_names, 2, mkChar("feval"));
    SET_STRING_ELT(out_names, 3, mkChar("intpval"));

    SET_NAMES(out, out_names);

    UNPROTECT(7);

    return out;
}

/* Added 2 by D Bates to move output from Powell Fortran codes */
 
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

void attribute_hidden F77_NAME(minqer)(const int *msgno)
{
    char *msg;
    switch(*msgno) {
    case 10:
	msg = _("bobyqa: NPT is not in the required interval");
	break;
    case 101:
	msg = _("newuoa: NPT is not in the required interval");
	break;
    case 20:
	msg = _("bobyqa: one of the differences XU(I)-XL(I) is less than 2*RHOBEG");
	break;
    case 320:
	msg = _("bobyqa detected too much cancellation in denominator");
	break;
    case 390:
	msg = _("bobyqa exceeded maximum number of function evaluations");
	break;
    case 430:
	msg = _("a trust region step in bobyqa failed to reduce q");
	break;
    case 3701:
	msg = _("a trust region step in newuoa failed to reduce q");
	break;
    case 2101:
	msg = _("a trust region step in uobyqa failed to reduce q");
	break;
    default:
	error(_("Unknown message number %d in f77err"), *msgno);
    }
    error(msg);
}

void attribute_hidden
F77_NAME(minqit)(const int *iprint, const double *rho, const int *nf,
		 const double *fopt, const int *n, const double xbase[],
		 const double xopt[])
{
    int i, ip = *iprint, nn = *n, nnf = *nf;
    if (ip >= 2) {
	Rprintf(_("Rho = %g, # of func. evals = %d\n"), *rho, *nf);
	Rprintf(_("Current min. f = %g at x of\n"), *fopt);
	for(i = 0; i < nn; i++) Rprintf("%g ", xbase[i] + xopt[i]);
	Rprintf("\n");
    }
}

void attribute_hidden
F77_NAME(minqi3)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint == 3) {
	int i, nn = *n;
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

void attribute_hidden
F77_NAME(minqir)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	int i, nn = *n;
	Rprintf(_("At return from bobyqa\n"));
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}
/* End 2 D Bates */

