#include "Rinternals.h"        // for SEXP

// for R_registerRoutines and R_CallMethodDef
#include "R_ext/Rdynload.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

  SEXP cgminu_wrapper(
      SEXP Rfun,
      SEXP Rgradient,
      SEXP Rinitial_values,
      SEXP Rcontrol);
  
  static R_CallMethodDef Rcgmin2_arg_description[] = {
    CALLDEF(cgminu_wrapper, 4),
    {NULL, NULL, 0} 
  };

  void R_init_Rcgmin2(DllInfo *info) {
    R_registerRoutines(info, NULL, Rcgmin2_arg_description, NULL, NULL);  
    R_useDynamicSymbols(info, FALSE);
  }

}  // extern "C"
