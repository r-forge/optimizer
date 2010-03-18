#include <Rcpp.h>
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>

/// Wrapper for the objective function.  It is initialized to R's "c" function
static Rcpp::Function cf("c");

/// Apply is.finite to elements in a transformation
static inline double check_finite(double x) {
    if (!R_finite(x)) Rf_error("non-finite parameter to user function");
    return x;
}

/// Fortran callable objective function evaluation.
extern "C"
double F77_NAME(calfun)(const int *n, const double x[], const int *ip) {
    Rcpp::Environment rho(cf.environment());
    Rcpp::IntegerVector cc(rho.get(".feval."));
    Rcpp::NumericVector pp(rho.get(".par."));

    cc[0]++;			// increment func eval count
    if (*n != pp.size())
	Rf_error("In calfun: n = %d but length(.par.) = %d", *n, pp.size());
    std::transform(x, x + *n, pp.begin(), check_finite); //also copies

    double f = Rcpp::as<double>(cf(pp)); // evaluate objective
    if (*ip == 3) {		// print eval info when very verbose
	Rprintf("%3d:%#14.8g:", cc[0], f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
    return f;
}

/// Names of returned values
//static Rcpp::Argument parnm("par"), fval("fval"), feval("feval");

/// Return the number of function evaluations as an SEXP
static SEXP fevalf() {
    Rcpp::Environment rho(cf.environment());
    return rho.get(".feval.");
}

/// Construct the classed list to return
static SEXP rval(Rcpp::NumericVector par, std::string cnm) {
    Rcpp::StringVector parnm(3), cl(2);
    int ip = 0, n = par.size();
    parnm = "par", "fval", "feval";
    cl = cnm, "minqa";

    Rcpp::List rr(3);
    rr = par, F77_NAME(calfun)(&n, par.begin(), &ip), fevalf();
    rr.names() = parnm;
    rr.attr("class") = cl;
    return wrap(rr);
}

/// Declaration of Powell's bobyqa
extern "C" 
void F77_NAME(bobyqa)(const int *n, const int *npt, double X[],
		      const double xl[], const double xu[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[]);

/// Interface for bobyqa
RcppExport SEXP bobyqa_cpp(SEXP ppar, SEXP pxl, SEXP pxu, SEXP pcc, SEXP fn) {
    Rcpp::Environment cc(pcc);
    Rcpp::NumericVector par(ppar), xl(pxl), xu(pxu);
    cf = Rcpp::Function(fn);	// install the objective function
    double rb = Rcpp::as<double>(cc.get("rhobeg")),
	re = Rcpp::as<double>(cc.get("rhoend"));
    int ip = Rcpp::as<int>(cc.get("iprint")),
	mxf = Rcpp::as<int>(cc.get("maxfun")),
	n = par.size(), np = Rcpp::as<int>(cc.get("npt"));
    std::vector<double> w((np + 5) * (np + n) + (3 * n * (n + 5))/2);
    
    F77_NAME(bobyqa)(&n, &np, par.begin(), xl.begin(), xu.begin(),
		     &rb, &re, &ip, &mxf, &w[0]);
    return rval(par, "bobyqa");
}

extern "C" 
void F77_NAME(uobyqa)(const int *n, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[]);

RcppExport SEXP uobyqa_cpp(SEXP ppar, SEXP pctrl, SEXP pfn) {
    Rcpp::Environment cc(pctrl); 
    Rcpp::NumericVector par(ppar);
    cf = Rcpp::Function(pfn);
    double rb = Rcpp::as<double>(cc.get("rhobeg")),
	re = Rcpp::as<double>(cc.get("rhoend"));
    int ip = Rcpp::as<int>(cc.get("iprint")),
	mxf = Rcpp::as<int>(cc.get("maxfun")), n = par.size();
    Rcpp::Environment rho(cf.environment());
    std::vector<double>
	w((n*(42+n*(23+n*(8+n))) + std::max(2*n*n + 4, 18*n)) / 4);

    F77_NAME(uobyqa)(&n, par.begin(), &rb, &re, &ip, &mxf, &w[0]);
    return rval(par, "uobyqa");
}

extern "C" 
void F77_NAME(newuoa)(const int *n, const int *npt, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[]);

RcppExport SEXP newuoa_cpp(SEXP ppar, SEXP pctrl, SEXP pfn) {
    Rcpp::Environment cc(pctrl); 
    Rcpp::NumericVector par(ppar);
    double rb = Rcpp::as<double>(cc.get("rhobeg")),
	re = Rcpp::as<double>(cc.get("rhoend"));
    int ip = Rcpp::as<int>(cc.get("iprint")),
	mxf = Rcpp::as<int>(cc.get("maxfun")),
	n = par.size(), np = Rcpp::as<int>(cc.get("npt"));
    std::vector<double> w((np+13)*(np+n)+(3*n*(n+3))/2);
    cf = Rcpp::Function(pfn);

    F77_NAME(newuoa)(&n, &np, par.begin(), &rb, &re, &ip, &mxf, &w[0]);
    return rval(par, "newuoa");
}

/// Assorted error messages.
extern "C" void F77_NAME(minqer)(const int *msgno) {
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
		 const double xopt[]) {
    if (*iprint >= 2) {
	Rprintf("%#8.2g; %3d: %#14.8g,", *rho, *nf, *fopt);
	for(int i = 0; i < *n; i++) Rprintf("%#8g ", xbase[i] + xopt[i]);
	Rprintf("\n");
    }
}

/// Output at return (do we really need this - why not use the print method?)
extern "C" void
F77_NAME(minqir)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[]) {
    if (*iprint > 0) {
	Rprintf("At return\n");
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(bobyqa_cpp, 5),
    CALLDEF(uobyqa_cpp, 3),
    CALLDEF(newuoa_cpp, 3),
    {NULL, NULL, 0}
};

/// Initializer for the package.  Registers the symbols for .Call.
extern "C" void R_init_minqa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



