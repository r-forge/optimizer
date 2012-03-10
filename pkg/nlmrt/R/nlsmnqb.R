nlsmnqb <-function(formula, start, trace=FALSE, data=NULL, 
            lower=-Inf, upper=Inf, masked=NULL, control=list(), ...){
#
#  A simplified and hopefully robust alternative to finding the 
#  nonlinear least squares minimizer that causes 'formula' to 
#  given a minimal residual sum of squares. 
#
# This version with bounds constraints
# ?? do we want masks e.g., lower>upper? Or explicitly??
#
#  nls.mn is particularly intended to allow for the resolution of 
#  very ill-conditioned or else near zero-residual problems for 
#  which the regular nls() function is ill-suited. It may also be
#  a useful pre-processor for nls().
#  
#  J C Nash  2012-3-4   nashjc _at_  uottawa.ca
#
#  formula looks like "y~b1/(1+b2*exp(-b3*T))"
#  start MUST be a vector where all the elements are named:
#   e.g., start=c(b1=200, b2=50, b3=0.3)
#  trace -- TRUE for extra output (?? NOT implemented)
# ?? control is a list of control parameters. These are:
#     ...
#  
#  ... will need to contain data for other variables that appear in
#  the formula and are defined in a parent frame
#
# data
#   dots<-list(...)
#   if (length(dots)>0){
#     print(dots)
#     str(dots)
#     xx<-..1
#     yy<-..2 # But how to get names right!??
#   }
# get data from data frame if exists
if (! is.null(data)){
    for (dfn in names(data)) {
       cmd<-paste(dfn,"<-data$",dfn,"")
       eval(parse(text=cmd))
    }
}
if (trace) {
  cat("Extracted data for variables:")
  print(names(data))
} 
# bounds
npar<-length(start) # number of parameters
if (length(lower)==1) lower<-rep(lower,npar)
if (length(upper)==1) upper<-rep(upper,npar) 
# ?? more tests on bounds / also scaling?
if (length(lower)!=npar) stop("Wrong length: lower")
if (length(upper)!=npar) stop("Wrong length: upper")
if (any(start<lower) || any(start>upper)) {
   cat("lower:")
   print(lower)
   cat("upper:")
   print(upper)
   stop("Infeasible start")
}
# controls
   ctrl<-list(
     watch=FALSE, # monitor progress
     phi=1, # the phi parameter
     lamda=0.0001, # lamda (spelled lamda in JNWMS)
     offset=100, # to determine if paramters changed
     laminc=10,
     lamdec=4 # use decreased_lamda<-lamda*lamdec/laminc
   )
   epstol<-(.Machine$double.eps)*ctrl$offset
   ncontrol <- names(control)
   nctrl <- names(ctrl)
   for (onename in ncontrol) {
      if (!(onename %in% nctrl)) {
         if (trace) cat("control ",onename," is not in default set\n")
      }
      ctrl[onename]<-control[onename]
   }
   phiroot<-sqrt(ctrl$phi)
   lamda<-ctrl$lamda # simplify controls to simplify code
   offset<-ctrl$offset
   laminc<-ctrl$laminc
   lamdec<-ctrl$lamdec # save typing
#  First get all the variable names:
    vn <- all.vars(parse(text=formula))
# Then see which ones are parameters (get their positions in the set xx
    pnum<-start # may simplify later??
    pnames<-names(pnum)
    bdmsk<-rep(1,npar) # set all params free for now
    bdmsk<-rep(1,npar) # set all params free for now
    mindx<-which(pnames==masked)
    bdmsk[mindx]<-0 # fixed parameters
    cat("bdmsk:")
    print(bdmsk)
    tmp<-readline("cont.")
    if (trace) {
      parpos  <- match(pnames, vn)
      datvar<-vn[-parpos] # NOT the parameters
      for (dvn in datvar){
         cat("Data variable ",dvn,":")
         print(eval(parse(text=dvn)))
      }
    }
    if (is.character(formula)){
       es<-formula
    } else {
       tstr<-as.character(formula) # note ordering of terms!
       es<-paste(tstr[[2]],"~",tstr[[3]],'')
    }
# Now separate the sides
    parts<-strsplit(as.character(es), "~")[[1]]
    if (length(parts)!=2) stop("Model expression is incorrect!")
    lhs<-parts[1]
    rhs<-parts[2]
# And build the residual at the parameters
    resexp<-paste(rhs,"-",lhs, collapse=" ")
    for (i in 1:npar){ # put parameters in separate variables
       joe<-paste(pnames[[i]],"<-",pnum[[i]])
       eval(parse(text=joe))
    }
    gradexp<-deriv(parse(text=resexp), names(start)) # gradient expression
    if (is.null(data)){
      resbest<-eval(parse(text=resexp)) # initial residual
    } else {resbest<-with(data, eval(parse(text=resexp))) }
    ssbest<-crossprod(resbest)
    feval<-1
    pbest<-pnum
    feval<-1 # number of function evaluations
    jeval<-0 # number of Jacobian evaluations
    if (trace) {
      cat("START -- lamda:",lamda," SS = ",ssbest," at ")
      print(pnum)
    }
    ssquares<-.Machine$double.xmax # make it big
    newjac<-TRUE # set control for computing Jacobian
    eqcount<-0
    while (eqcount < npar) {
       bdmsk[which(pnum-lower<=epstol*(abs(lower)+epstol))]<- -3 
       bdmsk[which(pnum-upper>=epstol*(abs(lower)+epstol))]<- -1
       cat("bdmsk:")
       print(bdmsk)
       if (newjac) {
          # cat("Evaluate Jacobian")
          if (is.null(data)){
            J0<-eval(gradexp) # initial residual
          } else {J0<-with(data, eval(gradexp))}
          Jac<-attr(J0,"gradient")
          jeval<-jeval+1 # count Jacobians
          if (any(is.na(Jac))) stop("NaN in Jacobian")
          dee<-diag(sqrt(diag(crossprod(Jac))))+phiroot*diag(npar) # to append to Jacobian
       }
       lamroot<-sqrt(lamda)
       JJ<-rbind(Jac,lamroot*dee) # build the matrix
       print(JJ)
       print(resbest)
       tmp<-readline("continue")
       JQR<-qr(JJ)# ??try
       rplus<-c(resbest, rep(0,npar))
       w<-qr.qty(JQR,rplus)
       gradproj<- -crossprod(w)
       cat("gradient projection=",gradproj,"\n") # ?? test
       delta<-try(qr.coef(JQR,-rplus)) # Note appended rows of y)
       if (class(delta)=="try-error") {
          if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
          lamda<-laminc*lamda
          newjac<-FALSE # increasing lamda -- don't re-evaluate
          if (trace) cat(" Singular matrix\n")
       } else {
#          if (trace) {
#             cat("delta:")
#             print(delta)
#          }
	  delta[mindx]<-0
          delta<-as.numeric(delta)
          cat("delta:")
          print(delta)
            step<-rep(1,npar)
            for (i in 1:npar){
                step[i]<-0
                if (delta[i]>0) step[i]<-(upper[i]-pbest[i])/delta[i]
                if (delta[i]<0) step[i]<-(lower[i]-pbest[i])/delta[i] # positive
            }
            print(step)
            stepsize<-min(1,step[which(delta!=0)])
            cat("Stepsize=",stepsize,"\n")
            tmp<-readline("more")
            pnum<-pbest+stepsize*delta # adjust (note POSITIVE here, but not in nlsmn0
            cat("pnum:")
            print(pnum)
#            names(pnum)<-pnames # NOT inherited through %*% !!!
	
          eqcount<-length(which((offset+pbest)==(offset+pnum)))
          if (eqcount<npar) {
             for (i in 1:npar){
               joe<-paste(pnames[[i]],"<-",pnum[[i]])
               eval(parse(text=joe))
             }
             feval<-feval+1 # count evaluations
             if (is.null(data)){
                resid<-eval(parse(text=resexp)) # trial residual
             } else {resid<-with(data, eval(parse(text=resexp)))}
             ssquares<-as.numeric(crossprod(resid))
             if (ssquares>=ssbest) {
                ## ?? see jnmws page 206 for simpler technique
                if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
                lamda<-laminc*lamda
                newjac<-FALSE # increasing lamda -- don't re-evaluate
                if(trace) cat(">= lamda=",lamda,"\n")
             } else {
                if (trace) {
                  cat("<< lamda=",lamda,"\n")
                  cat(" SS = ",ssquares," evals J/F:",jeval,"/",feval," eqcount=",eqcount,"\n")
                  print(pnum)
                }
                lamda<-lamdec*lamda/laminc
                ssbest<-ssquares
                resbest<-resid
                pbest<-pnum
                if (trace) {
                   cat("<< Lamda=",lamda,"\n")
                }
                newjac<-TRUE
             }
             if (ctrl$watch)  tmp<-readline() 
          } else {
             if (trace) cat(" No parameter change\n")
          }
      }
   }
## should return bdmsk too!
   result<-list(coeffs=pnum,ssquares=ssbest, resid=resbest, jacobian=Jac, feval=feval, jeval=jeval)
}


