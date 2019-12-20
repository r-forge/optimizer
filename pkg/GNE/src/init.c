/** 
 * @file  init.c
 * @brief init file for all C files
 *
 * @author Christophe Dutang
 *
 * Copyright (C) 2012, Christophe Dutang. 
 * All rights reserved.
 *  
 */
 
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "SSR.h"

 //table of registration for call methods
 static const R_CallMethodDef callMethods[] = 
   {
   {"dofunSSR", (DL_FUNC) &dofunSSR, 19},
   {"dojacSSR", (DL_FUNC) &dojacSSR, 24},
   {NULL, NULL, 0}
   };
 
 //table of registered routines
 void R_init_GNE(DllInfo *info)
 {
   //register method accessed with .Call
   R_registerRoutines(info, NULL, callMethods, NULL, NULL); 
   R_useDynamicSymbols(info, FALSE);
   R_forceSymbols(info, TRUE);
 }
 
