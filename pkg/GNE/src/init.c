/** 
 * @file  init.c
 * @brief init file for all C files
 *
 * @author Christophe Dutang
 *
 * Copyright (C) 2009, Christophe Dutang. 
 * All rights reserved.
 *  
 */
 
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "newton4KKT.h"
#include "util-mnl-pratio.h"
#include "util-mnl-pdiff.h"

//table of registration for call methods
static const R_CallMethodDef callMethods[] = 
{
        {"doPROBj2kRatio", (DL_FUNC) &doPROBj2kRatio, 4},
		{"doGradPROBj2kRatio", (DL_FUNC) &doGradPROBj2kRatio, 5},
		{"doGradGradPROBj2kRatio", (DL_FUNC) &doGradGradPROBj2kRatio, 6},
		{"doPROBj2kDiff", (DL_FUNC) &doPROBj2kDiff, 4},
		{"doGradPROBj2kDiff", (DL_FUNC) &doGradPROBj2kDiff, 5},
		{"doGradGradPROBj2kDiff", (DL_FUNC) &doGradGradPROBj2kDiff, 6},
		{"doPhi", (DL_FUNC) &doPhi, 12},
        {"doJacPhi", (DL_FUNC) &doJacPhi, 15},
        {NULL, NULL, 0}
};


//table of registered routines
void R_init_GNE(DllInfo *info)
{
        //register method accessed with .Call
        R_registerRoutines(info, NULL, callMethods, NULL, NULL); 
		
	
}

