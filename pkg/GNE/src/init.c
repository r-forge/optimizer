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
#include "util-clg.h"

//table of registration for call methods
static const R_CallMethodDef callMethods[] = 
{
        {"doPROBj2k", (DL_FUNC) &doPROBj2k, 4},
		{"doGradPROBj2k", (DL_FUNC) &doGradPROBj2k, 5},
		{"doGradGradPROBj2k", (DL_FUNC) &doGradGradPROBj2k, 6},
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

