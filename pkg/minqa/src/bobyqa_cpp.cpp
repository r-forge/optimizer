#include <Rcpp.h>
#include <R_ext/RS.h>

extern "C" void F77_NAME(minqer)(const int *msgno)
{
    const char *msg;
    switch(*msgno) {
    case 10:
	msg = "bobyqa: NPT is not in the required interval";
	break;
    case 101:
	msg = "newuoa: NPT is not in the required interval";
	break;
    case 20:
	msg = "bobyqa: one of the differences XU(I)-XL(I) is less than 2*RHOBEG";
	break;
    case 320:
	msg = "bobyqa detected too much cancellation in denominator";
	break;
    case 390:
	msg = "bobyqa exceeded maximum number of function evaluations";
	break;
    case 430:
	msg = "a trust region step in bobyqa failed to reduce q";
	break;
    case 3701:
	msg = "a trust region step in newuoa failed to reduce q";
	break;
    case 2101:
	msg = "a trust region step in uobyqa failed to reduce q";
	break;
    default:
	Rf_error("Unknown message number %d in f77err", *msgno);
    }
    Rf_error(msg);
}

extern "C" void
F77_NAME(minqit)(const int *iprint, const double *rho, const int *nf,
		 const double *fopt, const int *n, const double xbase[],
		 const double xopt[])
{
    int i, ip = *iprint, nn = *n;
    if (ip >= 2) {
	Rprintf("%#8.2g; %3d: %#14.8g,", *rho, *nf, *fopt);
	for(i = 0; i < nn; i++) Rprintf("%#8g ", xbase[i] + xopt[i]);
	Rprintf("\n");
    }
}

extern "C" void
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

extern "C" void
F77_NAME(minqir)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	int i, nn = *n;
	Rprintf("At return from bobyqa\n");
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

extern "C" void
F77_NAME(minqirn)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	int i, nn = *n;
	Rprintf("At return from newuoa\n");
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

extern "C" void
F77_NAME(minqiru)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[])
{
    if (*iprint > 0) {
	int i, nn = *n;
	Rprintf("At return from uobyqa\n");
	Rprintf("%3d:%#14.8g:", *nf, *f);
	for (i = 0; i < nn; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

Rcpp::Function cf("c");

extern "C" void F77_NAME(calfun)(const int *n, const double x[], double *f)
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

extern "C" 
void F77_NAME(bobyqa)(const int *n, const int *npt, double X[],
		      const double xl[], const double xu[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[],
		      double fval[]);

RcppExport SEXP bobyqa_cpp(SEXP par_arg, SEXP xl_arg, SEXP xu_arg, 
			   SEXP control, SEXP fn, SEXP work)
{
    Rcpp::NumericVector par(par_arg), xl(xl_arg), xu(xu_arg), wrk(work);
    Rcpp::Environment cc(control); 
    Rcpp::NumericVector rhobeg(cc.get("rhobeg")), rhoend(cc.get("rhoend"));
    Rcpp::IntegerVector npt(cc.get("npt")), maxfun(cc.get("maxfun")),
	iprint(cc.get("iprint"));
    int n = par.size();
    cf = Rcpp::Function(fn);

    if (!(rhobeg.size() && rhoend.size() && npt.size() && maxfun.size() &&
	  iprint.size())) Rf_error("Incomplete control parameter list");
    Rcpp::NumericVector fval(*(npt.begin()));
    Rcpp::NumericVector ff(1);

    if (xl.size() != n || xu.size() != n)
	Rf_error("sizes of par, xl and xu don't match");
    F77_NAME(bobyqa)(&n, npt.begin(), par.begin(), xl.begin(), xu.begin(),
		     rhobeg.begin(), rhoend.begin(), iprint.begin(),
		     maxfun.begin(), wrk.begin(), fval.begin());
    // For some bizarre reason the minimum function value is not returned
    F77_NAME(calfun)(&n, par.begin(), ff.begin());

    Rcpp::Environment rho(cf.environment());
    Rcpp::IntegerVector feval(rho.get(".feval."));
    Rcpp::List rr(Rcpp::Pairlist(Rcpp::Named("par", par),
				 Rcpp::Named("fval", ff),
				 Rcpp::Named("feval", feval),
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
