/** 
 * @file  newton4KKT.h
 * @brief header file for Newton type method to solve KKT system
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
SEXP doPhi(SEXP nplayer, SEXP z, SEXP dimx, SEXP dimlam, 
		   SEXP grobj, SEXP arggrobj, 
		   SEXP constr, SEXP argconstr, 
		   SEXP grconstr, SEXP arggrconstr, 
		   SEXP complfunc, SEXP env);


SEXP doJacPhi(SEXP nplayer, SEXP z, SEXP dimx, SEXP dimlam,
			  SEXP heobj, SEXP argheobj, 
			  SEXP constr, SEXP argconstr, 
			  SEXP grconstr, SEXP arggrconstr, 
			  SEXP heconstr, SEXP argheconstr,
			  SEXP grcompla, SEXP grcomplb, SEXP env);

 
