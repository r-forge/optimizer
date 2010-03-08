typedef struct opt_struct {
  SEXP par;
  SEXP fcall;
  SEXP env;
  double rss;
  int feval;
} opt_struct, *OptStruct;

SEXP getListElement(SEXP list, char *str);
double *real_vector(int n);
int  *int_vector(int n);

SEXP bobyqa_c(SEXP par_arg, SEXP xl, SEXP xu, SEXP fn, SEXP control, 
	      SEXP rho);
double resfunbobyqa(int n, double *par, double *f);

void F77_NAME(bobyqa)(int *N, int *NPT, double *X, double *XL, 
		      double *XU, double *RHOBEG, 
		      double *RHOEND, int *IPRINT, 
		      int *MAXFUN, double *W, double *FVAL);
		     
void F77_SUB(resfunbobyqa)(int *n, double *par, double *f); 
					   
extern OptStruct OS;


