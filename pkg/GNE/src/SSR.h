/** 
 * @file  newton4KKT.h
 * @brief header file for Newton type method to solve KKT system
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2012, Christophe Dutang.
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
SEXP dofunSSR(SEXP mode, SEXP nplayer, 
			  SEXP z, SEXP dimx, SEXP dimlam, SEXP dimmu,
			  SEXP grobj, SEXP arggrobj, 
			  SEXP constr, SEXP argconstr, 
			  SEXP grconstr, SEXP arggrconstr, 
			  SEXP joint, SEXP argjoint, 
			  SEXP grjoint, SEXP arggrjoint, 
			  SEXP complfunc, SEXP argcompl, SEXP env);


SEXP dojacSSR(SEXP mode, SEXP nplayer, 
			  SEXP z, SEXP dimx, SEXP dimlam, SEXP dimmu,
			  SEXP heobj, SEXP argheobj, 
			  SEXP constr, SEXP argconstr, 
			  SEXP grconstr, SEXP arggrconstr, 
			  SEXP heconstr, SEXP argheconstr,
			  SEXP joint, SEXP argjoint, 
			  SEXP grjoint, SEXP arggrjoint, 
			  SEXP hejoint, SEXP arghejoint,
			  SEXP grcompla, SEXP grcomplb, SEXP argcompl,
			  SEXP env);

 
