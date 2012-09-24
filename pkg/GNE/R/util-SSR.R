#functions of the SemiSmooth Reformulation of the GNEP
#z = (x, lambda, mu)


#function phi of the SSR
funSSR <- function(z, dimx, dimlam,
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	compl, argcompl, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	echo=FALSE)
{
#	cat("funSSR\n")
#	print(z)
#	print(dimx)
#	print(dimlam)
	
	arg <- testargfunSSR(z, dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
				  compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint, echo)

	mode <- ifelse(is.null(arg$constr), 0, 1) + 2*ifelse(is.null(arg$joint), 0, 1) 
#0 for no constr no joint
#1 for constr
#2 for joint
#3 for constr, joint
	
	res <- .Call("dofunSSR", as.integer(mode), arg$nplayer, z, 
				 as.integer(arg$dimx), as.integer(arg$dimlam), as.integer(arg$dimmu),
				 arg$grobj, arg$arggrobj, 
				 arg$constr, arg$argconstr, 
				 arg$grconstr, arg$arggrconstr, 
				 arg$joint, arg$argconstr, 
				 arg$grjoint, arg$arggrconstr, 
				 arg$compl, arg$argcompl, new.env())
	
	res
}


#Jacobian of phi of the SSR
jacSSR <- function(z, dimx, dimlam, 
	heobj, argheobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	hejoint, arghejoint,
	echo=FALSE)
{
	arg <- testargjacSSR(z, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
						 heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint,
						 hejoint, arghejoint, echo=FALSE)
	
	mode <- ifelse(is.null(arg$constr), 0, 1) + 2*ifelse(is.null(arg$joint), 0, 1) 
#0 for no constr no joint
#1 for constr
#2 for joint
#3 for constr, joint
	
	res <- .Call("dojacSSR", as.integer(mode), arg$nplayer, z, 
				 as.integer(arg$dimx), as.integer(arg$dimlam), as.integer(arg$dimmu),
				 arg$heobj, arg$argheobj, 
				 arg$constr, arg$argconstr, 
				 arg$grconstr, arg$arggrconstr, 
				 arg$heconstr, arg$argheconstr,
				 arg$joint, arg$argjoint, 
				 arg$grjoint, arg$arggrjoint, 
				 arg$hejoint, arg$arghejoint,
				 arg$gcompla, arg$gcomplb, arg$argcompl, new.env())

	
	res
}


