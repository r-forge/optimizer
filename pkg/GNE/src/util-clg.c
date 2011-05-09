/** 
 * @file  util-clg.c
 * @brief C file for cumulative logistic function
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2010, Christophe Dutang. 
 * All rights reserved.
 *
 */


#include "util-clg.h"




/*********************************/
/* cumulative logistic function */


//main function used .Call()
SEXP doPROBj2k(SEXP x, SEXP j, SEXP k, SEXP param_j)
{
    //temporary C working variables
	int i;
	int id_j, id_k; //indexes	
	int n; //nb element
	double *paramj; //parameters c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	double *c_x; //pointer to x
	double denom, num; //denominator and numerator
	double rmintilde; 

	//sanity check
    if (!isVector(x))
		error(_("invalid argument x"));
	if (length(x) <= 1)
		error(_("invalid length for argument x"));
	if (!isInteger(j))
		error(_("invalid argument j"));
	if (!isInteger(k))
		error(_("invalid argument k"));
	if (!isVector(param_j))
		error(_("invalid argument param_j"));	
	if (length(param_j) != 4)
		error(_("invalid length for argument param_j"));	
	
	//init
	id_j = asInteger(j) - 1; //index j
	id_k = asInteger(k) - 1; //index k	
	n = length(x); //nb element
	//parameters c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	paramj = REAL(param_j); 
	//pointer to x
	c_x = REAL(x); 

	
    //allocate result
    double *res; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocVector(REALSXP, 1)); //allocate a 1-vector
    res = REAL( resultinR ); //plug the C pointer on the R type
    
	
	//computation
	denom = 1.0;
	
	for (i = 0; i < n; i++)
		if (i != id_j && c_x[i] != 0.0)
			denom += exp(paramj[2] + paramj[3] * c_x[id_j] / c_x[i]);

	
	if (id_j == id_k)
	{
		num = 1.0;
		rmintilde = 1 - paramj[1];
	}else
	{
		num = exp(paramj[2] + paramj[3] * c_x[id_j] / c_x[id_k]);	
		rmintilde = paramj[0] / (n-1);
	}
	
	*res = rmintilde + (paramj[1] - paramj[0]) * num / denom;
	
	UNPROTECT(1);
	
	return resultinR;
}






