/** 
 * @file  util-clg.h
 * @brief header file for (util) multinomial logit (price-difference function)
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2010, Christophe Dutang.
 * All rights reserved.
 *
 */


//R header files
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Error.h>


#include "locale.h"



/* Functions accessed from .Call() */
SEXP doPROBj2kDiff(SEXP x, SEXP j, SEXP k, SEXP param_j);
SEXP doGradPROBj2kDiff(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP param_j);
SEXP doGradGradPROBj2kDiff(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP mderiv, SEXP param_j);

/* utility functions */
double PROBj2kDiff(double* x, int j, int k, double* param_j, int n);
double GradPROBj2kDiff(double* x, int j, int k, int ideriv, double* paramj, int n);
double GradGradPROBj2kDiff(double* x, int j, int k, int i, int m, double* paramj, int n); 
