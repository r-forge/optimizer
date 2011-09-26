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
/********************************/

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
    
	/*
	//computation - old code
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
	*/
	*res = PROBj2k(c_x, id_j, id_k, paramj, n);
	
	UNPROTECT(1);
	
	return resultinR;
}

//computation function by doPROBj2k
double PROBj2k(double* x, int j, int k, double* paramj, int n)
{
	//x price vector
	//j, k player indexes
	//paramj array of 4 elements c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	//n number of players
	
    //temporary C working variables
	int i;
	double denom, num; //denominator and numerator
	double rmintilde;  
	
	//computation
	denom = 1.0;
	
	for (i = 0; i < n; i++)
		if (i != j && x[i] != 0.0)
			denom += exp(paramj[2] + paramj[3] * x[j] / x[i]);
	
	if (j == k)
	{
		num = 1.0;
		rmintilde = 1 - paramj[1];
	}else
	{
		num = exp(paramj[2] + paramj[3] * x[j] / x[k]);	
		rmintilde = paramj[0] / (n-1);
	}

	return rmintilde + (paramj[1] - paramj[0]) * num / denom;
}



/*****************************************/
/* cumulative logistic gradient function */
/*****************************************/

//(x, j, k, ideriv, paramj)

//main function used .Call()
SEXP doGradPROBj2k(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP param_j)
{
    //temporary C working variables
	int id_j, id_k, id_i; //indexes	
	int n; //nb element
	double *paramj; //parameters c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	double *c_x; //pointer to x
	double denom, num; //denominator and numerator
	double rmintilde; 
	double *probarray; //array with prob(x, j, l) for l=1, ... n
	
	//sanity check
    if (!isVector(x))
		error(_("invalid argument x"));
	if (length(x) <= 1)
		error(_("invalid length for argument x"));
	if (!isInteger(j))
		error(_("invalid argument j"));
	if (!isInteger(k))
		error(_("invalid argument k"));
	if (!isInteger(ideriv))
		error(_("invalid argument ideriv"));
	if (!isVector(param_j))
		error(_("invalid argument param_j"));	
	if (length(param_j) != 4)
		error(_("invalid length for argument param_j"));	
	
	//init
	id_j = asInteger(j) - 1; //index j
	id_k = asInteger(k) - 1; //index k	
	id_i = asInteger(ideriv) - 1; //derivative index i
	n = length(x); //nb element
	//parameters e.g. c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	paramj = REAL(param_j); 
	//pointer to x
	c_x = REAL(x); 
	
	
    //allocate result
    double *res; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocVector(REALSXP, 1)); //allocate a 1-vector
    res = REAL( resultinR ); //plug the C pointer on the R type
    
		//(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP param_j)
	*res = GradPROBj2k(c_x, id_j, id_k, id_i, paramj, n);
	
	UNPROTECT(1);
	
	return resultinR;
}

//computation function by doGradPROBj2k
double GradPROBj2k(double* x, int j, int k, int i, double* paramj, int n)
{
	//x price vector
	//j, k player indexes
	//i derivative index
	//paramj array of 4 elements c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	//n number of players
	
	//temporary C working variables
	int l;
	double Sumlclgxjl; //array of sum_l clg(x, j, l) / x_l for l=1, ..., n	and l diff j
	double clgxjk; //clg(x, j, k) 
	
	//result 
	double res;
	
	clgxjk = PROBj2k(x, j, k, paramj, n);
	
	if(i == j)
	{	
		
		Sumlclgxjl = 0.0;
		for (l = 0; l < n; l++)
			if(l != j && x[l] != 0.0)
				Sumlclgxjl += PROBj2k(x, j, l, paramj, n) / x[l];
		
		res = -paramj[3] * Sumlclgxjl * clgxjk;
	}else
	{
		if(x[i] != 0.0)
			res = paramj[3] * x[j] / (x[i] * x[i]) * PROBj2k(x, j, i, paramj, n) * clgxjk; 
		else
			res = 0.0;
	}
	
	if(i == j && j != k)
		res += paramj[3] / x[k] * clgxjk;  
	
	if(i == k && j != k && x[k] != 0.0)
		res -= paramj[3] * x[j] / (x[k] * x[k]) * clgxjk; 
	
	return res;
}	



/*****************************************/
/* cumulative logistic Hessian function */
/*****************************************/


//(x, j, k, ideriv, mderiv, paramj)

//main function used .Call()
SEXP doGradGradPROBj2k(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP mderiv, SEXP param_j)
{
    //temporary C working variables
	int id_j, id_k, id_i, id_m; //indexes	
	int n; //nb element
	double *paramj; //parameters c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	double *c_x; //pointer to x
	double denom, num; //denominator and numerator
	double rmintilde; 
	double *probarray; //array with prob(x, j, l) for l=1, ... n
	
	//sanity check
    if (!isVector(x))
		error(_("invalid argument x"));
	if (length(x) <= 1)
		error(_("invalid length for argument x"));
	if (!isInteger(j))
		error(_("invalid argument j"));
	if (!isInteger(k))
		error(_("invalid argument k"));
	if (!isInteger(ideriv))
		error(_("invalid argument ideriv"));
	if (!isInteger(mderiv))
		error(_("invalid argument mderiv"));
	if (!isVector(param_j))
		error(_("invalid argument param_j"));	
	if (length(param_j) != 4)
		error(_("invalid length for argument param_j"));	
	
	//init
	id_j = asInteger(j) - 1; //index j
	id_k = asInteger(k) - 1; //index k	
	id_i = asInteger(ideriv) - 1; //derivative index i
	id_m = asInteger(mderiv) - 1; //derivative index i	
	n = length(x); //nb element
	//parameters e.g. c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	paramj = REAL(param_j); 
	//pointer to x
	c_x = REAL(x); 
	
	
    //allocate result
    double *res; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocVector(REALSXP, 1)); //allocate a 1-vector
    res = REAL( resultinR ); //plug the C pointer on the R type
    
	//(SEXP x, SEXP j, SEXP k, SEXP ideriv, SEXP param_j)
	*res = GradGradPROBj2k(c_x, id_j, id_k, id_i, id_m, paramj, n);
	
	UNPROTECT(1);
	
	return resultinR;
}

//computation function by doGradGradPROBj2k
double GradGradPROBj2k(double* x, int j, int k, int i, int m, double* paramj, int n)
{
	//x price vector
	//j, k player indexes
	//i, m derivative indexes
	//paramj array of 4 elements c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
	//n number of players
	
	//temporary C working variables
	int l;
	double Sumlclgxjl; //array of sum_l clg(x, j, l) / x_l for l=1, ..., n	and l diff j
	double Sumlgrmclgxjl; //array of sum_l Gr x_m clg(x, j, l) / x_l for l=1, ..., n	and l diff j	
	double clgxjk, grmclgxjk;
	double clgxji;
	
	double res = 0.0;
	
	clgxjk = PROBj2k(x, j, k, paramj, n);
	grmclgxjk = GradPROBj2k(x, j, k, m, paramj, n);
	
	
	if(i == j)
	{
		Sumlclgxjl = 0.0;
		for (l = 0; l < n; l++)
			if(l != j && x[l] != 0.0)
				Sumlclgxjl += PROBj2k(x, j, l, paramj, n) / x[l];
		
		Sumlgrmclgxjl = 0.0;
		for (l = 0; l < n; l++)
			if(l != j && x[l] != 0.0)
				Sumlgrmclgxjl += GradPROBj2k(x, j, l, m, paramj, n) / x[l];
		
		if(x[k] != 0.0)
			res += -paramj[3] * Sumlgrmclgxjl * clgxjk;
		res += -paramj[3] * Sumlclgxjl * grmclgxjk;
		if(j != m && x[m] != 0.0)
			res += paramj[3] / (x[m] * x[m]) * clgxjk * PROBj2k(x, j, m, paramj, n);
	}else
	{
		clgxji = PROBj2k(x, j, i, paramj, n);
		
		if(x[i] != 0.0)
		{	
			res += paramj[3] * x[j] / (x[i] * x[i]) * GradPROBj2k(x, j, i, m, paramj, n) * clgxjk;
			res += paramj[3] * x[j] / (x[i] * x[i]) * clgxji * grmclgxjk;
			if(j == m)
				res += paramj[3] / (x[i] * x[i]) * clgxji * clgxjk;
			if(i == m)
				res += -2*paramj[3] * x[j] / (x[i] * x[i] * x[i]) * clgxji * clgxjk; 
		}	
	}
	
	
	if(j != k && i == j)
	{
		if(x[k] != 0.0)
		{
			res += paramj[3] / x[k] * grmclgxjk;
			if(k == m)	
				res += -paramj[3] / (x[m] * x[m]) * clgxjk;
		}
	}
	
	if(j !=k && i == k && x[i] != 0.0)
	{
		res += -paramj[3] * x[j] / (x[i] * x[i]) * grmclgxjk;
		if(j == m)
			res += -paramj[3] / (x[i] * x[i]) * clgxjk;
		if(i == m)
			res += 2*paramj[3] * x[j] / (x[i] * x[i] * x[i]) * clgxjk;
	}
	
	return res;
}



