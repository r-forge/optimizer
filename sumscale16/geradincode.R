<<geradincode, echo=FALSE, cache=TRUE>>=
ax<-function(x, AA){
   u<-as.numeric(AA%*%x)
}

bx<-function(x, BB){
   v<-as.numeric(BB%*%x)
}

geradin<-function(x, ax, bx, AA, BB, control=list(trace=TRUE, maxit=1000)){
# Geradin minimize Rayleigh Quotient, Nash CMN Alg 25
  print(control)
  trace<-control$trace
  n<-length(x)
  tol<-n*n*.Machine$double.eps^2
  offset<-1e+5 # equality check offset
  if (trace) cat("geradin.R, using tol=",tol,"\n")
  ipr<-0 # counter for matrix mults
  pa<-.Machine$double.xmax
  R<-pa
  msg<-"no msg"
# step 1 -- main loop
  keepgoing<-TRUE
  while (keepgoing) {
    avec<-ax(x, AA); bvec<-bx(x, BB); ipr<-ipr+1
    xax<-as.numeric(crossprod(x, avec));  
    xbx<-as.numeric(crossprod(x, bvec));
    if (xbx <= tol) {
       keepgoing<-FALSE # not really needed
       msg<-"avoid division by 0 as xbx too small"
       break
    } 
    p0<-xax/xbx
    if (p0>pa) {
       keepgoing<-FALSE # not really needed
       msg<-"Rayleigh Quotient increased in step"
       break
    } 
    pa<-p0
    g<-2*(avec-p0*bvec)/xbx
    gg<-as.numeric(crossprod(g)) # step 6
    if (trace) cat("Before loop: RQ=",p0," after ",ipr," products, gg=",gg,"\n")
    if (gg<tol) { # step 7
       keepgoing<-FALSE # not really needed
       msg<-"Small gradient -- done"
       break
    } 
    t<- -g # step 8
    for (itn in 1:n) { # major loop step 9
       y<-ax(t, AA); z<-bx(t, BB); ipr<-ipr+1 # step 10
       tat<-as.numeric(crossprod(t, y)) # step 11
       xat<-as.numeric(crossprod(x, y)) 
       xbt<-as.numeric(crossprod(x, z)) 
       tbt<-as.numeric(crossprod(t, z)) 
       u<-tat*xbt-xat*tbt
       v<-tat*xbx-xax*tbt
       w<-xat*xbx-xax*xbt
       d<-v*v-4*u*w
       if (d<0) stop("Geradin: imaginary roots not possible") # step 13
       d<-sqrt(d) # step 14
       if (v>0) k<--2*w/(v+d) else k<-0.5*(d-v)/u
       xlast<-x # NOT as in CNM -- can be avoided with loop
       avec<-avec+k*y; bvec<-bvec+k*z # step 15, update
       x<-x+k*t
       xax<-xax+as.numeric(crossprod(x,avec))      
       xbx<-xbx+as.numeric(crossprod(x,bvec))      
       if (xbx<tol) stop("Geradin: xbx has become too small")
       chcount<-n - length(which((xlast+offset)==(x+offset)))
       if (trace) cat("Number of changed components = ",chcount,"\n")
       pn<-xax/xbx # step 17 different order
       if (chcount==0) {
         keepgoing<-FALSE # not really needed
         msg<-"Unchanged parameters -- done"
         break
       }
       if (pn >= p0) {
         if (trace) cat("RQ not reduced, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       p0<-pn # step 19
       g<-2*(avec-pn*bvec)/xbx
       gg<-as.numeric(crossprod(g))
       if (trace) cat("Itn", itn," RQ=",p0," after ",ipr," products, gg=",gg,"\n")
       if (gg<tol){ # step 20
         if (trace) cat("Small gradient in iteration, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       xbt<-as.numeric(crossprod(x,z)) # step 21
       w<-y-pn*z # step 22
       tabt<-as.numeric(crossprod(t,w))
       beta<-as.numeric(crossprod(g,(w-xbt*g)))
       beta<-beta/tabt # step 23
       t<-beta*t-g
    } # end loop on itn -- step 24
  } # end main loop -- step 25
# step 26
  ans<-list(x=x, RQ=p0, ipr=ipr, msg=msg)
}
@

<<rungeradin1, echo=FALSE, cache=TRUE>>=

cat("Test geradin\n")
n<-10
AA<-molerbuild(n)
x<-runif(n)
#cat("A Matrix:\n")
#print(AA)
BB=diag(rep(1,n))
#cat("B Matrix:\n")
#print(BB)
cat("xstart:")
print(x)
tg<-system.time(ag<-geradin(x, ax, bx, AA=AA, BB=BB, 
   control=list(trace=FALSE)))[[3]]
cat("Minimal eigensolution\n")
ag$x<-ag$x/sqrt(as.numeric(crossprod(ag$x))) # rescale
print(ag)
cat("Geradin time=",tg,"\n")
tgn<-system.time(agn<-geradin(x, ax, bx, AA=-AA, BB=BB,
   control=list(trace=FALSE)))[[3]]
cat("Maximal eigensolution (negative matrix)\n")
agn$x<-agn$x/sqrt(as.numeric(crossprod(agn$x))) # rescale
print(agn)
cat("Geradin time=",tgn,"\n")
@

