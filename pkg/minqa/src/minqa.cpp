#include <Rcpp.h>
#include <R_ext/RS.h>

/// Wrapper for the objective function.  It is initialized to R's "c" function
Rcpp::Function cf("c");

/// Fortran callable objective subroutine.  Why is this not a function?
extern "C"
void F77_NAME(calfun)(const int *n, const double x[], double *f)
{
    Rcpp::Environment rho(cf.environment());
    Rcpp::IntegerVector cc(rho.get(".feval."));
    Rcpp::NumericVector pp(rho.get(".par."));
    double *p = pp.begin();

    (*cc.begin())++;		// increment the function evaluation count
    if (*n != pp.size())
	Rf_error("In calfun: n = %d but length(.par.) = %d", *n, pp.size());
    for (int i = 0; i < *n; i++) {
	if (!R_finite(p[i]))
	    Rf_error("non-finite parameter to user function");
	p[i] = x[i];
    }
    *f = Rcpp::NumericVector(cf()).begin()[0];
}

/// Declaration of Powell's bobyqa
extern "C" 
void F77_NAME(bobyqa)(const int *n, const int *npt, double X[],
		      const double xl[], const double xu[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[],
		      double fval[]);

/// Interface for bobyqa
RcppExport SEXP bobyqa_cpp(SEXP ppar, SEXP pxl, SEXP pxu, 
			   SEXP ctrl, SEXP fn, SEXP work)
{
    Rcpp::Environment cc(ctrl); 
    Rcpp::NumericVector par(ppar), xl(pxl), xu(pxu), wrk(work),
	rb(cc.get("rhobeg")), re(cc.get("rhoend"));
    int n = par.size();
    if (xl.size() != n || xu.size() != n) // check sizes (lengths)
	Rf_error("sizes of par, xl and xu don't match");
    Rcpp::IntegerVector npt(cc.get("npt")), mxf(cc.get("maxfun")),
	ip(cc.get("iprint"));

    cf = Rcpp::Function(fn);	// install the objective function
    Rcpp::Environment rho(cf.environment());

    // Ensure that all the control settings contain at least one element
    if (!(rb.size() && re.size() && npt.size() && mxf.size() && ip.size()))
	Rf_error("Incomplete control parameter list");

    // Allocate memory for function values used for approximation.  Is
    // this used at all? It looks like it was added after the fact.
    Rcpp::NumericVector fval(*(npt.begin())); 

    F77_NAME(bobyqa)(&n, npt.begin(), par.begin(), xl.begin(),
		     xu.begin(), rb.begin(), re.begin(), ip.begin(),
		     mxf.begin(), wrk.begin(), fval.begin());

    // For some bizarre reason the minimum function value is not returned
    Rcpp::NumericVector ff(1);
    F77_NAME(calfun)(&n, par.begin(), ff.begin());

    Rcpp::List rr(Rcpp::Pairlist(Rcpp::Named("par", par),
				 Rcpp::Named("fval", ff),
				 Rcpp::Named("feval", rho.get(".feval.")),
				 Rcpp::Named("intpval", fval)));
    rr.attr("class") = "bobyqa";
    return rr;
}

extern "C" 
void F77_NAME(uobyqa)(const int *n, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[]);

RcppExport SEXP uobyqa_cpp(SEXP ppar, SEXP pctrl, SEXP pfn, SEXP work)
{
    Rcpp::NumericVector par(ppar), wrk(work), ff(1);
    Rcpp::Environment cc(pctrl); 
    Rcpp::NumericVector rhobeg(cc.get("rhobeg")), rhoend(cc.get("rhoend"));
    Rcpp::IntegerVector maxfun(cc.get("maxfun")), iprint(cc.get("iprint"));
    cf = Rcpp::Function(pfn);
    int n = par.size();

    if (!(rhobeg.size() && rhoend.size() && maxfun.size() && iprint.size()))
	Rf_error("Incomplete control parameter list");

    F77_NAME(uobyqa)(&n, par.begin(), rhobeg.begin(), rhoend.begin(),
		     iprint.begin(), maxfun.begin(), wrk.begin());
    // For some bizarre reason the minimum function value is not returned
    F77_NAME(calfun)(&n, par.begin(), ff.begin());

    Rcpp::Environment rho(cf.environment());
    Rcpp::IntegerVector feval(rho.get(".feval."));
    Rcpp::List rr(Rcpp::Pairlist(Rcpp::Named("par", par),
				 Rcpp::Named("fval", ff),
				 Rcpp::Named("feval", feval)));
    rr.attr("class") = "uobyqa";
    return rr;
}

extern "C" 
void F77_NAME(newuoa)(const int *n, const int *npt, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[]);

RcppExport SEXP newuoa_cpp(SEXP ppar, SEXP pctrl, SEXP pfn, SEXP work)
{
    Rcpp::NumericVector par(ppar), wrk(work), ff(1);
    Rcpp::Environment cc(pctrl); 
    Rcpp::NumericVector rhobeg(cc.get("rhobeg")), rhoend(cc.get("rhoend"));
    Rcpp::IntegerVector maxfun(cc.get("maxfun")), iprint(cc.get("iprint")),
	npt(cc.get("npt"));
    cf = Rcpp::Function(pfn);

    if (!(rhobeg.size() && npt.size() && rhoend.size() &&
	  maxfun.size() && iprint.size()))
	Rf_error("Incomplete control parameter list");

    int n = par.size();
    F77_NAME(newuoa)(&n, npt.begin(), par.begin(), rhobeg.begin(),
		     rhoend.begin(), iprint.begin(), maxfun.begin(),
		     wrk.begin());
    // For some bizarre reason the minimum function value is not returned
    F77_NAME(calfun)(&n, par.begin(), ff.begin());

    Rcpp::Environment rho(cf.environment());
    Rcpp::IntegerVector feval(rho.get(".feval."));
    Rcpp::List rr(Rcpp::Pairlist(Rcpp::Named("par", par),
				 Rcpp::Named("fval", ff),
				 Rcpp::Named("feval", feval)));
    rr.attr("class") = "uobyqa";
    return rr;
}

/// Assorted error messages.
extern "C" void F77_NAME(minqer)(const int *msgno)
{
    const char *msg;
    switch(*msgno) {
    case 10:
    case 101:
	msg = "NPT is not in the required interval";
	break;
    case 20:
	msg = "one of the differences XU(I)-XL(I) is less than 2*RHOBEG";
	break;
    case 320:
	msg = "bobyqa detected too much cancellation in denominator";
	break;
    case 390:
	msg = "maximum number of function evaluations exceeded";
	break;
    case 430:
    case 3701:
    case 2101:
	msg = "a trust region step failed to reduce q";
	break;
    default:
	Rf_error("Unknown message number %d in minqer", *msgno);
    }
    Rf_error(msg);
}

/// Iteration output when rho changes and iprint >= 2
extern "C" void
F77_NAME(minqit)(const int *iprint, const double *rho, const int *nf,
		 const double *fopt, const int *n, const double xbase[],
		 const double xopt[])
{
    if (*iprint >= 2) {
	Rprintf("%#8.2g; %3d: %#14.8g,", *rho, *nf, *fopt);
	for(int i = 0; i < *n; i++) Rprintf("%#8g ", xbase[i] + xopt[i]);
	Rprintf("\n");
    }
}

/// Function evaluation trace output for iprint == 3
extern "C" void
F77_NAME(minqi3)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint == 3) {
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

/// Output at return (do we really need this - why not use the print method?)
extern "C" void
F77_NAME(minqir)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	Rprintf("At return\n");
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}