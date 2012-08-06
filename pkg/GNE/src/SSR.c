/** 
 * @file  funSSR.c
 * @brief C file for functions of the SemiSmooth Reformulation of the GNEP
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2012, Christophe Dutang. 
 * All rights reserved.
 *
 */


#include "SSR.h"




/*********************************/
/* function phi of the SSR */


//main function used .Call()
SEXP dofunSSR(SEXP mode, SEXP nplayer, 
			  SEXP z, SEXP dimx, SEXP dimlam, SEXP dimmu,
			  SEXP grobj, SEXP arggrobj, 
			  SEXP constr, SEXP argconstr, 
			  SEXP grconstr, SEXP arggrconstr, 
			  SEXP joint, SEXP argjoint, 
			  SEXP grjoint, SEXP arggrjoint, 
			  SEXP complfunc, SEXP argcompl, SEXP env)
{
    if (!isInteger(mode))
        error(_("invalid argument mode"));
    if (!isVector(z))
        error(_("invalid argument z"));
	if (!isFunction(grobj))
		error(_("invalid argument grobj"));
	if (!isFunction(complfunc))
		error(_("invalid argument complfunc"));
    
    //temporary C working variables
	int MODE = asInteger(mode);
	int i, j, k, idz, idlam; // index;
    int N = asInteger(nplayer); //number of players
	int n, m, max_m; //dimension of x, dimension of lambda, maximum dimension of lambda_i
	int cumsum_j;	//cumsum_i
	int *dim_x = INTEGER(dimx);
	int *dim_lam = INTEGER(dimlam);
	int dim_mu = asInteger(dimmu);
	double *lambda, *mu;
	
	if (MODE == 1 || MODE == 3)
		if (!isFunction(constr) || !isFunction(grconstr))
			error(_("invalid argument constr, grconstr"));	
	if (MODE == 2 || MODE == 3)
		if (!isFunction(joint) || !isFunction(grjoint))
			error(_("invalid argument joint, grjoint"));	
	
	
    //declare temporary R working variables
	SEXP play; //variable for the player number
	SEXP deriv1; //index for the 1st derivative
	SEXP R_grobjcall, grobj_res; //gradient of the objective function
	SEXP R_constrcall, constr_res; //the constraint function
	SEXP R_grconstrcall, grconstr_res; //gradient of the constraint function
	SEXP R_jointcall, joint_res; //the joint function
	SEXP R_grjointcall, grjoint_res; //gradient of the joint function
	SEXP R_complcall, compl_res; //the complementarity function
	SEXP a,b; //arguments of compl func
	
	//init C variables
	n = 0; m = 0; max_m = 0;
	for(i = 0; i < N; i++) 
		n += dim_x[i];
	if (MODE == 1 || MODE == 3)
	{
		for(i = 0; i < N; i++) 
		{
			m += dim_lam[i];
			if(dim_lam[i] > max_m)
				max_m = dim_lam[i];
		}
		
		lambda  = (double *) R_alloc(m, sizeof(double));
		for(i = 0; i < m; i++) 
			lambda[i] = REAL(z)[n+i];
	}
	if (MODE == 2 || MODE == 3)
	{
		mu  = (double *) R_alloc(dim_mu, sizeof(double));
		for(i = 0; i < dim_mu; i++) 
			mu[i] = REAL(z)[n+m+i];
	}
	
	idz = 0; //index in z
	idlam = 0; //index in lambda
	//for idz <= n, z[idz] = x[idz] and 
	//for m+n > idz > n, z[idz] = lambda[idz-n]
	//for idz > n+m, z[idz] = mu[idz-n-m]
	
	if(length(z) != n && MODE == 0)
		error(_("problem of dimension for z."));		
	if(length(z) != n+m && MODE == 1)
		error(_("problem of dimension for z."));		
	if(length(z) != n+dim_mu && MODE == 2)
		error(_("problem of dimension for z."));		
	if(length(z) != n+m+dim_mu && MODE == 3)
		error(_("problem of dimension for z."));		
	
//		Rprintf("check n:%d m:%d l:%d \n", n, m, dim_mu);
	
	//alloc SEXP variables for function
	PROTECT(play = allocVector(INTSXP, 1));
	PROTECT(deriv1 = allocVector(INTSXP, 1));		
	PROTECT(a = allocVector(REALSXP, 1));
	PROTECT(b = allocVector(REALSXP, 1));
	
	//alloc SEXP variables for function calls
	PROTECT(R_grobjcall = lang5(grobj, z, play, deriv1, arggrobj)); //Grad_i O_i(x) 
	PROTECT(R_complcall = lang4(complfunc, a, b, argcompl)); //phi(a,b)   
	if (MODE == 1 || MODE == 3)
	{
		PROTECT(R_constrcall = lang4(constr, z, play, argconstr)); //g_i(x) 
		PROTECT(R_grconstrcall = lang5(grconstr, z, play, deriv1, arggrconstr)); //Grad_j g_i(x) 	
	}
	if (MODE == 2 || MODE == 3)
	{
		PROTECT(R_jointcall = lang3(joint, z, argjoint)); //h(x) 
		PROTECT(R_grjointcall = lang4(grjoint, z, deriv1, arggrjoint)); //Grad_j h(x) 	
	}
	
	//alloc SEXP variables for results
	PROTECT(grobj_res = allocVector(REALSXP, 1)); 
	PROTECT(compl_res = allocVector(REALSXP, 1));
	if (MODE == 1 || MODE == 3)
	{
		PROTECT(constr_res = allocVector(REALSXP, max_m));
		PROTECT(grconstr_res = allocVector(REALSXP, max_m));	
	}
	if (MODE == 2 || MODE == 3)
	{
		PROTECT(joint_res = allocVector(REALSXP, dim_mu));
		PROTECT(grjoint_res = allocVector(REALSXP, dim_mu));	
	}
	
	
    //alloc result
    double *phiz; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocVector(REALSXP, n+m+dim_mu)); //allocate a (n+m+dim_mu) vector
    phiz = REAL( resultinR ); //plug the C pointer on the R type

    R_CheckStack();
	
	/**************************************************
	 * computation
	 *
	 * Phi = (A)
	 *       (B) 
	 *       (C) 	 
	 **************************************************/

	/***************************/
	//computation of the A-part

	//computation of Grad_i L_i = Grad_i O_i(x) + sum_k lambda[i]^T_k Grad_i (g_i(x))_k + mu^T Grad_i (h(x))
	cumsum_j = 0;
	for(i = 0; i < N; i++) //player i
	{	
		INTEGER(play)[0] = i+1;
		
		for(j = 0; j < dim_x[i]; j++)	//index for variable x	
		{
			INTEGER(deriv1)[0] = cumsum_j+j+1;
			PROTECT(grobj_res = eval(R_grobjcall, env)); //Grad_{x_ij}  O_i(x) 
			phiz[idz] = REAL(grobj_res)[0];

			if (MODE == 1 || MODE == 3)
			{
				PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_ij} g_i(x)
				
				for(k = 0; k < dim_lam[i]; k++)	//index for variable lambda	
					phiz[idz] += lambda[idlam + k] * REAL(grconstr_res)[k];
				UNPROTECT(1);
			}
			if (MODE == 2 || MODE == 3)
			{
				PROTECT(grjoint_res = eval(R_grjointcall, env)); //Grad_{x_ij} h(x)
				
				for(k = 0; k < dim_mu; k++)	//index for variable mu	
					phiz[idz] += mu[k] * REAL(grjoint_res)[k];
				UNPROTECT(1);
			}
			idz++;	
			UNPROTECT(1);
		}
		cumsum_j += dim_x[i];
		idlam += dim_lam[i];		
	}	
	
	/***************************/
	//computation of the B-part
 	
	//computation of phi_i = phi(-g_i(x), lambda[i]) component wise
	if (MODE == 1 || MODE == 3)
	{
		for(i = 0; i < N; i++)
		{	
			INTEGER(play)[0] = i+1;
			PROTECT(constr_res = eval(R_constrcall, env)); //g_i(x)

			for(j = 0; j < dim_lam[i]; j++)
			{
//				Rprintf("test 2 -REAL(constr_res)[%d]:%f \n", j, -REAL(constr_res)[j]);
				REAL(a)[0] = -REAL(constr_res)[j];
				REAL(b)[0] = REAL(z)[idz]; //lambda[idz - n]
				PROTECT(compl_res = eval(R_complcall, env)); //phi(-g_i(x), lambda[i])
				phiz[idz] = REAL(compl_res)[0];				
				idz++;	 
			}
			UNPROTECT(1+dim_lam[i]);
		}
	}
	
		
	
	/***************************/
	//computation of the C-part
 	
	//computation of phi(-h(x), mu) component wise
	if (MODE == 2 || MODE == 3)
	{
		PROTECT(joint_res = eval(R_jointcall, env)); //h(x)
			
		for(j = 0; j < dim_mu; j++)
		{
//			Rprintf("test 3 -REAL(joint_res)[%d]:%f \n", j, -REAL(joint_res)[j]);			
			REAL(a)[0] = -REAL(joint_res)[j];
			REAL(b)[0] = REAL(z)[idz]; //mu[idz - n - m]
			PROTECT(compl_res = eval(R_complcall, env)); //phi(-h(x), mu)
			phiz[idz] = REAL(compl_res)[0];				
			idz++;	 
		}
		UNPROTECT(1+dim_mu);
	}
	
    UNPROTECT(9);
	if (MODE == 1 || MODE == 2)
		UNPROTECT(4);
	if (MODE == 3)
		UNPROTECT(8);

    return resultinR;
}




/*********************************/
/* Jacobian of phi of the SSR */

//main function used .Call()
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
			  SEXP env)
{
	if (!isInteger(mode))
        error(_("invalid argument mode"));
    if (!isVector(z))
        error(_("invalid argument z"));
	if (!isFunction(heobj))
		error(_("invalid argument heobj"));
	if (!isFunction(grcompla) || !isFunction(grcomplb))
		error(_("invalid argument grcompla, grcomplb"));
    
    
    //temporary C working variables
	int MODE = asInteger(mode);
	int i, j, k, l, idlam; // idz, idgrconstr;
	int cumsum_i, cumsum_j;
    int N = asInteger(nplayer); //number of players
	int n, m, max_m; //dimension of x, dimension of lambda, maximum dimension of lambda_i
	int *dim_x = INTEGER(dimx);
	int *dim_lam = INTEGER(dimlam);
	int dim_mu = asInteger(dimmu); 
	double *lambda, *mu;
	double *tempconstr, *tempjoint;
	
	if (MODE == 1 || MODE == 3)
		if (!isFunction(constr) || !isFunction(grconstr))
			error(_("invalid argument constr, grconstr"));	
	if (MODE == 2 || MODE == 3)
		if (!isFunction(joint) || !isFunction(grjoint))
			error(_("invalid argument joint, grjoint"));	
	
    //init C variables
	n = 0; m = 0; max_m = 0;
	for(i = 0; i < N; i++) 
		n += dim_x[i];
	if (MODE == 1 || MODE == 3)
	{
		for(i = 0; i < N; i++) 
		{
			m += dim_lam[i];
			if(dim_lam[i] > max_m)
				max_m = dim_lam[i];
		}
		lambda = (double *) R_alloc(m, sizeof(double));
		
		for(i = 0; i < m; i++) 
			lambda[i] = REAL(z)[n+i];
		tempconstr = (double *) R_alloc(m, sizeof(double));
	}
	if (MODE == 2 || MODE == 3)
	{
		mu  = (double *) R_alloc(dim_mu, sizeof(double));
		for(i = 0; i < dim_mu; i++) 
			mu[i] = REAL(z)[n+m+i];
		tempjoint = (double *) R_alloc(dim_mu, sizeof(double));
	}
	R_CheckStack();
	
	if(length(z) != n && MODE == 0)
		error(_("problem of dimension for z."));		
	if(length(z) != n+m && MODE == 1)
		error(_("problem of dimension for z."));		
	if(length(z) != n+dim_mu && MODE == 2)
		error(_("problem of dimension for z."));		
	if(length(z) != n+m+dim_mu && MODE == 3)
		error(_("problem of dimension for z."));	
	
	//declare R working variables
	SEXP play; //player variable
	SEXP deriv1; //index for first derivative
	SEXP deriv2; //index for first derivative
	SEXP R_heobjcall, heobj_res; //the Hessian of the obj func
	SEXP R_constrcall, constr_res; //the constraint function
	SEXP R_grconstrcall, grconstr_res; //gradient of constr func
	SEXP R_heconstrcall, heconstr_res; //the Hessian of constr func
	SEXP R_jointcall, joint_res; //the joint function
	SEXP R_grjointcall, grjoint_res; //gradient of joint func
	SEXP R_hejointcall, hejoint_res; //the Hessian of joint func
	SEXP R_grcomplacall, grcompla_res; //first derivative of compl func phi'_a
	SEXP R_grcomplbcall, grcomplb_res; //first derivative of compl func phi'_b	
	SEXP a,b; //arguments of compl func a,b
	
	
	//alloc SEXP variables for function
	PROTECT(play = allocVector(INTSXP, 1));	
	PROTECT(deriv1 = allocVector(INTSXP, 1));	
	PROTECT(deriv2 = allocVector(INTSXP, 1));	
	PROTECT(a = allocVector(REALSXP, 1));
	PROTECT(b = allocVector(REALSXP, 1));
	
	//alloc SEXP variables for function calls
	PROTECT(R_heobjcall = lang6(heobj, z, play, deriv1, deriv2, argheobj)); //Grad_k Grad_j 0_i(x)
	PROTECT(R_grcomplacall = lang4(grcompla, a, b, argcompl)); //phi'_a(a,b)
	PROTECT(R_grcomplbcall = lang4(grcomplb, a, b, argcompl)); //phi'_b(a,b)	
	if (MODE == 1 || MODE == 3)
	{
		PROTECT(R_constrcall = lang4(constr, z, play, argconstr)); //g_i(x) 
		PROTECT(R_grconstrcall = lang5(grconstr, z, play, deriv1, arggrconstr)); //Grad_j g_i(x) 	
		PROTECT(R_heconstrcall = lang6(heconstr, z, play, deriv1, deriv2, argheconstr)); //Grad_k Grad_j g_i(x)
	}
	if (MODE == 2 || MODE == 3)
	{
		PROTECT(R_jointcall = lang3(joint, z, argjoint)); //h(x) 
		PROTECT(R_grjointcall = lang4(grjoint, z, deriv1, arggrjoint)); //Grad_j h(x) 	
		PROTECT(R_hejointcall = lang5(hejoint, z, deriv1, deriv2, arghejoint)); //Grad_k Grad_j h(x)
	}
	
	//alloc SEXP variables for results
	PROTECT(heobj_res = allocVector(REALSXP, 1)); 
	PROTECT(grcompla_res = allocVector(REALSXP, 1));
	PROTECT(grcomplb_res = allocVector(REALSXP, 1));	
	if (MODE == 1 || MODE == 3)
	{
		PROTECT(constr_res = allocVector(REALSXP, max_m)); 
		PROTECT(grconstr_res = allocVector(REALSXP, max_m));
		PROTECT(heconstr_res = allocVector(REALSXP, max_m));	
	}
	if (MODE == 2 || MODE == 3)
	{
		PROTECT(joint_res = allocVector(REALSXP, dim_mu)); 
		PROTECT(grjoint_res = allocVector(REALSXP, dim_mu));
		PROTECT(hejoint_res = allocVector(REALSXP, dim_mu));	
	}
	
    //alloc result
    double *jacphiz; //result in C
    SEXP resultinR; //result in R
    PROTECT(resultinR = allocMatrix(REALSXP, n+m+dim_mu, n+m+dim_mu)); //alloc a (n+m++dim_mu) x (n+m++dim_mu) matrix
    jacphiz = REAL( resultinR ); //plug the C pointer on the R type
	//beware jacphiz is stored column by column
	
	/**************************************************
	 * init
	 **************************************************/    
	
    //init the result matrix 
	for(i = 0; i < n+m+dim_mu; i++)
		for(j = 0; j < n+m+dim_mu; j++)
			jacphiz[i + j * (n+m+dim_mu)] = 0.0;
	
/*	for(i = 0; i < n+m+dim_mu; i++) Rprintf("z: %f \n", REAL(z)[i]); Rprintf("\n");
	for(i = 0; i < m; i++) Rprintf("lam: %f \n", lambda[i]); Rprintf("\n");	
	for(i = 0; i < dim_mu; i++) Rprintf("mu: %f \n", mu[i]); Rprintf("\n");	
*/
	
	if (MODE == 1 || MODE == 3)
	{
		idlam = 0;
		//compute the constr function g_1(x) ... g_N(x)
		for(i = 0; i < N; i++)
		{	
			INTEGER(play)[0] = i+1;
			PROTECT(constr_res = eval(R_constrcall, env)); //g_i(x)
			for(j = 0; j < dim_lam[i]; j++)
				tempconstr[idlam + j] = REAL(constr_res)[j];
			idlam += dim_lam[i];		
			UNPROTECT(1);
		}
		/* for(i = 0; i < m; i++) Rprintf("constr: %f \n", tempconstr[i]); Rprintf("\n"); */	
	}

	if (MODE == 2 || MODE == 3)
	{
		//compute the joint function h(x)
		PROTECT(joint_res = eval(R_jointcall, env)); //h(x)
		for(j = 0; j < dim_mu; j++)
			tempjoint[j] = REAL(joint_res)[j];
		UNPROTECT(1);
		/* for(i = 0; i < dim_mu; i++) Rprintf("joint: %f \n", tempjoint[i]); Rprintf("\n"); */
	}	
	
	/**************************************************
	 * computation
	 *
	 * JacPhi = (A B C)
	 *          (D E F)
	 *			(G H I)
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
				//	Rprintf("player %d - deriv index i %d, j %d\t\t", l+1, i+1, j+1);
				jacphiz[i + j*(n+m+dim_mu)] = REAL(heobj_res)[0]; 				
				UNPROTECT(1);

				if (MODE == 1 || MODE == 3)
				{
					PROTECT(heconstr_res = eval(R_heconstrcall, env)); //Grad_{x_j} Grad_{x_i} g_l(x)
					for(k = 0; k < dim_lam[l]; k++)
					{	
						jacphiz[i + j*(n+m+dim_mu)] += lambda[idlam + k] * REAL(heconstr_res)[k];
						//Rprintf("k %d lambda %f he constr %f\t", k, lambda[idlam + k], REAL(heconstr_res)[k]);
					}
					UNPROTECT(1);
				}
				if (MODE == 2 || MODE == 3)
				{
					PROTECT(hejoint_res = eval(R_hejointcall, env)); //Grad_{x_j} Grad_{x_i} h(x)
					for(k = 0; k < dim_mu; k++)
					{
						jacphiz[i + j*(n+m+dim_mu)] += mu[k] * REAL(hejoint_res)[k];
						//Rprintf("k %d mu %f he joint %f\t", k+1, mu[k], REAL(hejoint_res)[k]);
					}
					UNPROTECT(1);
				}
			}
		}
		cumsum_i += dim_x[l];
		idlam += dim_lam[l];		
	}	
	
	/***************************/
	//computation of the B-part	
	if (MODE == 1 || MODE == 3)
	{
		cumsum_j = n; //cumsum index for j's
		cumsum_i = 0; //cumsum index for i's
		for(l = 0; l < N; l++) //player
		{	
			INTEGER(play)[0] = l+1;
			//		Rprintf("player: %d \n cum sum i %d \t cum sum j %d\n", l+1, cumsum_i, cumsum_j);
			for(i = cumsum_i; i < cumsum_i + dim_x[l]; i++) //row index
			{	
				INTEGER(deriv1)[0] = i+1;		
				PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_li} g_l(x)
				
				for(j = cumsum_j; j < cumsum_j + dim_lam[l]; j++) //column index
				{
					jacphiz[i + j*(n+m+dim_mu)] = REAL(grconstr_res)[j-cumsum_j];
				}
			}
			UNPROTECT(dim_x[l]);
			cumsum_i += dim_x[l];
			cumsum_j += dim_lam[l];
		}
	}

	/***************************/
	//computation of the C-part	
	if (MODE == 2 || MODE == 3)
	{
		cumsum_j = n+m; //cumsum index for j's
		cumsum_i = 0; //cumsum index for i's
		for(l = 0; l < N; l++) //player
		{
			for(i = cumsum_i; i < cumsum_i + dim_x[l]; i++) //row index
			{	
				INTEGER(deriv1)[0] = i+1;		
				PROTECT(grjoint_res = eval(R_grjointcall, env)); //Grad_{x_li} h(x)
				
				for(j = cumsum_j; j < cumsum_j + dim_mu; j++) //column index
				{
					jacphiz[i + j*(n+m+dim_mu)] = REAL(grjoint_res)[j-cumsum_j];
				}
			}
			UNPROTECT(dim_x[l]);
			cumsum_i += dim_x[l];		
		}			
	}
		
	/***************************/
	//computation of the D-part	
	if (MODE == 1 || MODE == 3)
	{
		cumsum_i = n; //cumsum index for i's
		k = 0;
		for(l = 0; l < N; l++) //player
		{	
			INTEGER(play)[0] = l+1;
			for(j = 0; j < n; j++) //column index
			{	
				INTEGER(deriv1)[0] = j+1;		
				PROTECT(grconstr_res = eval(R_grconstrcall, env)); //Grad_{x_lj} g_l(x)
				for(i = cumsum_i; i < cumsum_i + dim_lam[l]; i++) //row index
				{
					k++;
					REAL(a)[0] = -tempconstr[i-n]; //-g_i(x)
					REAL(b)[0] = REAL(z)[i]; //lambda[i]
					PROTECT(grcompla_res = eval(R_grcomplacall, env)); //phi'_a(-g_i(x), lambda[i])			
					//	 Rprintf("i %d, j %d, phi'a %f\n", i, j, REAL(grcompla_res)[0]);
					jacphiz[i + j*(n+m+dim_mu)] = - REAL(grconstr_res)[i-cumsum_i] * REAL(grcompla_res)[0];
				}
				UNPROTECT(1+dim_lam[l]);
			}
			cumsum_i += dim_lam[l];
		}
	}
	
	
	/***************************/
	//computation of the E-part
	if (MODE == 1 || MODE == 3)
	{
		for(i = n; i < n+m; i++)
		{	
			REAL(a)[0] = -tempconstr[i-n]; //-g_i(x)
			REAL(b)[0] = REAL(z)[i]; //lambda[i]
			PROTECT(grcomplb_res = eval(R_grcomplbcall, env)); //phi'_b(-g_i(x), lambda[i])			
			jacphiz[i + i*(n+m+dim_mu)] = REAL(grcomplb_res)[0];	
			UNPROTECT(1);
		}
	}
	
	/***************************/
	//computation of the F-part
	/* already to zero */	
	
	/***************************/
	//computation of the G-part	
	if (MODE == 2 || MODE == 3)
	{
		cumsum_i = n+m; //cumsum index for i's
		k = 0;
		for(j = 0; j < n; j++) //column index
		{	
			INTEGER(deriv1)[0] = j+1;		
			PROTECT(grjoint_res = eval(R_grjointcall, env)); //Grad_{x_lj} h(x)
			//			Rprintf("index derivative %d\n", j+1);
			for(i = cumsum_i; i < cumsum_i + dim_mu; i++) //row index
			{
				k++;
				REAL(a)[0] = -tempjoint[i-cumsum_i]; //-h(x)
				REAL(b)[0] = REAL(z)[i]; //mu[i]
				PROTECT(grcompla_res = eval(R_grcomplacall, env)); //phi'_a(-h(x), mu[i])			
				jacphiz[i + j*(n+m+dim_mu)] = - REAL(grjoint_res)[i-cumsum_i] * REAL(grcompla_res)[0];
				//				Rprintf("i %d, j %d, phi'a %f\n", i, j, REAL(grcompla_res)[0]);					
			}
			UNPROTECT(1+dim_mu);
		}
	}	
	
	/***************************/
	//computation of the H-part
	/* already to zero */
	
	/***************************/
	//computation of the I-part
	if (MODE == 2 || MODE == 3)
	{
		for(i = n+m; i < n+m+dim_mu; i++)
		{	
			REAL(a)[0] = -tempjoint[i-n-m]; //-h(x)
			REAL(b)[0] = REAL(z)[i]; //mu[i]
			PROTECT(grcomplb_res = eval(R_grcomplbcall, env)); //phi'_b(-h(x), mu[i])			
			jacphiz[i + i*(n+m+dim_mu)] = REAL(grcomplb_res)[0];	
			UNPROTECT(1);
		}
	}
	
	
	UNPROTECT(12);
	if (MODE == 1 || MODE == 2)
		UNPROTECT(6);
	if (MODE == 3)
		UNPROTECT(12);
	
    return resultinR;
}







