/** 
 * @file  newton4KKT.c
 * @brief C file for Newton type methods to solve KKT system
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2010, Christophe Dutang. 
 * All rights reserved.
 *
 */


#include "newton4KKT.h"




/*********************************/
/* Phi function */


//main function used .Call()
SEXP doPhi(SEXP nplayer, SEXP z, SEXP dimx, SEXP dimlam, 
		   SEXP grobj, SEXP arggrobj, 
		   SEXP constr, SEXP argconstr, 
		   SEXP grconstr, SEXP arggrconstr, 
		   SEXP complfunc, SEXP env)
{
//    if (!isReal(nplayer))
  //      error(_("invalid argument nplayer"));
    if (!isVector(z))
        error(_("invalid argument z"));
	if (!isFunction(grobj) || !isFunction(constr) || !isFunction(grconstr))
		error(_("invalid argument grobj, constr, grconstr"));
	if (!isFunction(complfunc))
		error(_("invalid argument complfunc"));
	
//	if (!isList(arggrobj))
  //      error(_("invalid argument arggrobj"));

    
    //temporary C working variables
	int i, j, k, idz, idlam; // idgrconstr;
    int N = asInteger(nplayer); //number of players
	int n, m; //dimension of x and dimension of lambda
	int cumsum_j;	//cumsum_i
	int *dim_x = INTEGER(dimx);
	int *dim_lam = INTEGER(dimlam);
//	double *grobj_res_C, *constr_res_C;
	double *lambda;

	
    //temporary R working variables
	SEXP play; //variable for the player number
	SEXP deriv1; //index for the 1st derivative
	SEXP R_grobjcall, grobj_res; //gradient of the objective function
	SEXP R_constrcall, constr_res; //the constraint function
	SEXP R_grconstrcall, grconstr_res; //gradient of the constraint function
	SEXP R_complcall, compl_res; //the complementarity function
	SEXP a,b; //arguments of compl func
	
	//init C variables
	n = 0; m = 0; 
	for(i = 0; i < N; i++) 
		n += dim_x[i];
	for(i = 0; i < N; i++) 
		m += dim_lam[i];

	lambda  = (double *) R_alloc(m, sizeof(double));
	for(i = 0; i < m; i++) 
		lambda[i] = REAL(z)[n+i];

	idz = 0; //index in z
	//for idz <= n, z[idz] = x[idz] and for idz > n, z[idz] = lambda[idz-n]
	idlam = 0; //index in lambda

	
	if(length(z) != n+m)
		error(_("problem of dimension"));		
	
	
	//init SEXP variables for function
	PROTECT(play = allocVector(INTSXP, 1));
	PROTECT(deriv1 = allocVector(INTSXP, 1));		
	PROTECT(a = allocVector(REALSXP, 1));
	PROTECT(b = allocVector(REALSXP, 1));
	
	//init SEXP variables for function calls
	PROTECT(R_grobjcall = lang5(grobj, z, play, deriv1, arggrobj)); //Grad_i O_i(x) 
	PROTECT(R_constrcall = lang4(constr, z, play, argconstr)); //g_i(x) 
	PROTECT(R_grconstrcall = lang5(grconstr, z, play, deriv1, arggrconstr)); //Grad_j g_i(x) 	
	PROTECT(R_complcall = lang3(complfunc, a, b)); //phi(a,b)   
	
	//init SEXP variables for results
	PROTECT(grobj_res = allocVector(REALSXP, 1)); 
	PROTECT(constr_res = allocVector(REALSXP, m));
	PROTECT(grconstr_res = allocVector(REALSXP, m));	
	PROTECT(compl_res = allocVector(REALSXP, 1));
	
    //allocate result
    double *phiz; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocVector(REALSXP, n+m)); //allocate a (n+m) vector
    phiz = REAL( resultinR ); //plug the C pointer on the R type

    R_CheckStack();


	
	/**************************************************
	 * computation
	 *
	 * Phi = (A)
	 *       (B) 
	 **************************************************/
	

	/***************************/
	//computation of the A-part

	//computation of Grad_i L_i = Grad_i O_i(x) + sum_k lambda[i]^T_k Grad_i (g_i(x))_k
	cumsum_j = 0;
	for(i = 0; i < N; i++) //player i
	{	
		INTEGER(play)[0] = i+1;
		
		for(j = 0; j < dim_x[i]; j++)	//index for variable x	
		{
			INTEGER(deriv1)[0] = cumsum_j+j+1;
			
			PROTECT(grobj_res = eval(R_grobjcall, env)); //Grad_{x_ij}  O_i(x) 
			
			PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_ij} g_i(x)
			
			phiz[idz] = REAL(grobj_res)[0];

//			Rprintf("res %f\n", phiz[idz] );			
//			Rprintf("player i %d - idx deriv %d - dim %d \n", i+1, cumsum_j+j+1, dim_x[i]);
			for(k = 0; k < dim_lam[i]; k++)	//index for variable lambda
			{

				phiz[idz] += lambda[idlam + k] * REAL(grconstr_res)[k];
//				Rprintf("i %d, j %d, k %d \n", i,j,k);
//				Rprintf("res %f\n", phiz[idz] );
//				Rprintf("lam_k: %f \n", lambda[idlam + k]);
//				Rprintf("gr_i g_i: %f \n", REAL(grconstr_res)[k]);				
			}
			
//			Rprintf("\n");
			idz++;	
			
			UNPROTECT(2);

		}
		cumsum_j += dim_x[i];
		idlam += dim_lam[i];		
	}
	

	/***************************/
	//computation of the B-part
 	
	//computation of phi_i = phi(-g_i(x), lambda[i]) component wise
	for(i = 0; i < N; i++)
	{	
		INTEGER(play)[0] = i+1;
		PROTECT(constr_res = eval(R_constrcall, env)); //g_i(x)

//				Rprintf("z %f\n", REAL(z)[idz] );			
		
		for(j = 0; j < dim_lam[i]; j++)
		{
			REAL(a)[0] = -REAL(constr_res)[j];
			REAL(b)[0] = REAL(z)[idz]; //lambda[idz - n]
			PROTECT(compl_res = eval(R_complcall, env)); //phi(-g_i(x), lambda[i])

//				Rprintf("compl %f\n", compl_res );						
			
			phiz[idz] = REAL(compl_res)[0];
			
			idz++;	 
		}
		UNPROTECT(1+dim_lam[i]);
	}
	
    UNPROTECT(13);

    return resultinR;
}




/*********************************/
/* Jacobian of the Phi function */



//main function used .Call()
SEXP doJacPhi(SEXP nplayer, SEXP z, SEXP dimx, SEXP dimlam,
			  SEXP heobj, SEXP argheobj, 
			  SEXP constr, SEXP argconstr, 
			  SEXP grconstr, SEXP arggrconstr, 
			  SEXP heconstr, SEXP argheconstr,
			  SEXP grcompla, SEXP grcomplb, SEXP env)
{
    if (!isVector(z))
        error(_("invalid argument z"));
	if (!isFunction(heobj) || !isFunction(constr) || !isFunction(grconstr) || !isFunction(heconstr))
		error(_("invalid argument heobj, constr, grconstr, heconstr"));
	if (!isFunction(grcompla) || !isFunction(grcomplb))
		error(_("invalid argument grcompla, grcomplb"));
    
    
    //temporary C working variables
	int i, j, k, l, idlam; // idz, idgrconstr;
	int cumsum_i, cumsum_j;
    int N = asInteger(nplayer); //number of players
	int n, m, max_m; //dimension of x, dimension of lambda, maximum dimension of lambda_i
	int *dim_x = INTEGER(dimx);
	int *dim_lam = INTEGER(dimlam);
	double *lambda;
	double *tempconstr;
	
    //init C variables
	n = 0; m = 0; 
	for(i = 0; i < N; i++) 
		n += dim_x[i];
	for(i = 0; i < N; i++) 
		m += dim_lam[i];
	max_m = dim_lam[0];
	for(i = 0; i < N; i++) 
		if(dim_lam[i] > max_m)
			max_m = dim_lam[i];
	lambda = (double *) R_alloc(m, sizeof(double));
	R_CheckStack();
	for(i = 0; i < m; i++) 
		lambda[i] = REAL(z)[n+i];
	tempconstr = (double *) R_alloc(m, sizeof(double));
	R_CheckStack();
	
	if(length(z) != n+m)
		error(_("problem of dimension"));		

	
	//temporary R working variables
	SEXP play; //player variable
	SEXP deriv1; //index for first derivative
	SEXP deriv2; //index for first derivative
	SEXP R_heobjcall, heobj_res; //the Hessian of the obj func
	SEXP R_constrcall, constr_res; //the constraint function
	SEXP R_grconstrcall, grconstr_res; //gradient of constr func
	SEXP R_heconstrcall, heconstr_res; //the Hessian of constr func
	SEXP R_grcomplacall, grcompla_res; //first derivative of compl func phi'_a
	SEXP R_grcomplbcall, grcomplb_res; //first derivative of compl func phi'_b	
	SEXP a,b; //arguments of compl func a,b
	
	
	//init SEXP variables for function
	PROTECT(play = allocVector(INTSXP, 1));	
	PROTECT(deriv1 = allocVector(INTSXP, 1));	
	PROTECT(deriv2 = allocVector(INTSXP, 1));	
	PROTECT(a = allocVector(REALSXP, 1));
	PROTECT(b = allocVector(REALSXP, 1));
	
	//init SEXP variables for function calls
	PROTECT(R_heobjcall = lang6(heobj, z, play, deriv1, deriv2, argheobj)); //Grad_k Grad_j 0_i(x)
	PROTECT(R_constrcall = lang4(constr, z, play, argconstr)); //g_i(x) 
	PROTECT(R_grconstrcall = lang5(grconstr, z, play, deriv1, arggrconstr)); //Grad_j g_i(x) 	
	PROTECT(R_heconstrcall = lang6(heconstr, z, play, deriv1, deriv2, argheobj)); //Grad_k Grad_j g_i(x)
	PROTECT(R_grcomplacall = lang3(grcompla, a, b)); //phi'_a(a,b)
	PROTECT(R_grcomplbcall = lang3(grcomplb, a, b)); //phi'_b(a,b)	
	
	//init SEXP variables for results
	PROTECT(heobj_res = allocVector(REALSXP, 1)); 
	PROTECT(constr_res = allocVector(REALSXP, max_m)); 
	PROTECT(grconstr_res = allocVector(REALSXP, max_m));
	PROTECT(heconstr_res = allocVector(REALSXP, max_m));	
	PROTECT(grcompla_res = allocVector(REALSXP, 1));
	PROTECT(grcomplb_res = allocVector(REALSXP, 1));	
	
    //allocate result
    double *jacphiz; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocMatrix(REALSXP, n+m, n+m)); //allocate a (n+m) x (n+m) matrix
    jacphiz = REAL( resultinR ); //plug the C pointer on the R type
	//beware jacphiz is stored column by column
	
    
	
    //init the result matrix 
	for(i = 0; i < n+m; i++)
		for(j = 0; j < n+m; j++)
			jacphiz[i + j * (n+m)] = 0.0;
	
/*	
	for(i = 0; i < n+m; i++)
		Rprintf("z: %f \n", REAL(z)[i]);
	Rprintf("\n");
	for(i = 0; i < m; i++)
		Rprintf("lam: %f \n", lambda[i]);
	Rprintf("\n");	
	
	Rprintf("max m: %d \n", max_m);
*/
	
	idlam = 0;
	//compute the constr function g_1(x) ... g_N(x)
	for(i = 0; i < N; i++)
	{	
		INTEGER(play)[0] = i+1;
//		Rprintf("play: %d\n", INTEGER(play)[0]);
		PROTECT(constr_res = eval(R_constrcall, env)); //g_i(x)
		for(j = 0; j < dim_lam[i]; j++)
		{	
//			Rprintf("j: %d , idlam+j: %d,\t\t g: %f\n", j, idlam + j, REAL(constr_res)[j]);
			tempconstr[idlam + j] = REAL(constr_res)[j];
		}
		idlam += dim_lam[i];		
		
		UNPROTECT(1);
	}
	
/*	
	for(i = 0; i < m; i++)
		Rprintf("constr: %f \n", tempconstr[i]);
	Rprintf("\n");
*/	
	
	/**************************************************
	 * computation
	 *
	 * JacPhi = (A B)
	 *          (C D) 
	 **************************************************/


	
	/***************************/
	//computation of the A-part
	cumsum_i = 0; //cumsum index for i's
	idlam = 0; //index for lambda
	for(l = 0; l < N; l++) //player
	{	
		INTEGER(play)[0] = l+1;
		
		for(i = cumsum_i; i < cumsum_i + dim_x[l]; i++) //row index
		{
			INTEGER(deriv1)[0] = i+1;		
			
			for(j = 0; j < n; j++) //column index
			{

				
				INTEGER(deriv2)[0] = j+1;		
				PROTECT(heobj_res = eval(R_heobjcall, env)); //Grad_{x_j} Grad_{x_i} O_l(x)
				
//				Rprintf("player %d - deriv index i %d, j %d\t\t", l+1, i+1, j+1);
				
				jacphiz[i + j*(n+m)] = REAL(heobj_res)[0]; 
				
				PROTECT(heconstr_res = eval(R_heconstrcall, env)); //Grad_{x_j} Grad_{x_i} g_l(x)
				
				for(k=0; k < dim_lam[l]; k++)
				{	
					jacphiz[i + j*(n+m)] += lambda[idlam + k] * REAL(heconstr_res)[k];
//					Rprintf("idx lambda %d, lambda %f\n", idlam+k, lambda[idlam + k]);
				}
				
				UNPROTECT(2);
			}
			
		}
		cumsum_i += dim_x[l];
		idlam += dim_lam[l];		
	}	
	
	/***************************/
	//computation of the B-part	
	cumsum_j = n; //cumsum index for j's
	cumsum_i = 0; //cumsum index for i's
	for(l = 0; l < N; l++) //player
	{	
		INTEGER(play)[0] = l+1;
		
//		Rprintf("player: %d \n", l+1);
//		Rprintf("cum sum i %d \t cum sum j %d\n", cumsum_i, cumsum_j);
		for(i = cumsum_i; i < cumsum_i + dim_x[l]; i++) //row index
		{	
			INTEGER(deriv1)[0] = i+1;		
			PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_li} g_l(x)
			
			for(j = cumsum_j; j < cumsum_j + dim_lam[l]; j++) //column index
			{
				jacphiz[i + j*(n+m)] = REAL(grconstr_res)[j-cumsum_j];
			}
		}
		cumsum_i += dim_x[l];
		cumsum_j += dim_lam[l];
	}
	
	k=0;
	
	/***************************/
	//computation of the C-part	
	cumsum_i = n; //cumsum index for i's
	
	for(l = 0; l < N; l++) //player
	{	
		INTEGER(play)[0] = l+1;
		
		for(j = 0; j < n; j++) //column index
		{	

			INTEGER(deriv1)[0] = j+1;		
			PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_lj} g_l(x)
//			Rprintf("index derivative %d\n", j+1);
			
			for(i = cumsum_i; i < cumsum_i + dim_lam[l]; i++) //row index
			{
				k++;
				REAL(a)[0] = -tempconstr[i-n]; //-g_i(x)
				REAL(b)[0] = REAL(z)[i]; //lambda[i]
//				Rprintf("a: %f, b: %f \n", -tempconstr[i-n], REAL(z)[i]);
				PROTECT(grcompla_res = eval(R_grcomplacall, env)); //phi'_a(-g_i(x), lambda[i])			
				
//				Rprintf("i %d, j %d, phi'a %f\n", i, j, REAL(grcompla_res)[0]);
				
				jacphiz[i + j*(n+m)] = - REAL(grconstr_res)[i-cumsum_i] * REAL(grcompla_res)[0];
				
//				Rprintf("k: %d \n", k);
			}
			UNPROTECT(1+dim_lam[l]);
		}
		cumsum_i += dim_lam[l];
	}
	
	
	/***************************/
	//computation of the D-part
	for(i = n; i < n+m; i++)
	{	
		REAL(a)[0] = -tempconstr[i-n]; //-g_i(x)
		REAL(b)[0] = REAL(z)[i]; //lambda[i]
		PROTECT(grcomplb_res = eval(R_grcomplbcall, env)); //phi'_b(-g_i(x), lambda[i])			
		jacphiz[i + i * (n+m)] = REAL(grcomplb_res)[0];	
	}
	
	
    UNPROTECT(18 + n + m);
	

	
    return resultinR;
}







