## ----setup, echo=FALSE, cache=FALSE--------------------------------------
# size takes valid value of LaTeX font sizes like small, big, huge, ...
# opts_chunk$set(size = 'scriptsize', fig.pos='htbp')
library(knitr)
opts_chunk$set(size = 'scriptsize')

## ----molermat, echo=TRUE-------------------------------------------------
molermat<-function(n){
   A<-matrix(NA, nrow=n, ncol=n)
   for (i in 1:n){
      for (j in 1:n) {
          if (i == j) A[i,i]<-i
          else A[i,j]<-min(i,j) - 2
      }
   }
   A
}

## ----molerfast, echo=TRUE------------------------------------------------
molerfast <- function(n) {
# A fast version of `molermat'
A <- matrix(0, nrow = n, ncol = n)
j <- 1:n
for (i in 1:n) {
A[i, 1:i] <- pmin(i, 1:i) - 2
}
A <- A + t(A)
diag(A) <- 1:n
A
}

## ----molertime, echo=FALSE-----------------------------------------------
nmax<-10
mtable<-matrix(NA, nrow=nmax, ncol=7) # to hold results
require(compiler)
molerc<-cmpfun(molermat) # compile it
molerfc<-cmpfun(molerfast)
# loop over sizes
for (ni in 1:nmax){
  n<-50*ni
  mtable[[ni, 1]]<-n
  ti<-system.time(ai<-molermat(n))[[1]]
  tc<-system.time(ac<-molerc(n))[[1]]
  if (! identical(ai, ac)) stop("Different outcomes == molermat, molerc")
  tfi<-system.time(afi<-molerfast(n))[[1]]
  tfc<-system.time(afc<-molerfc(n))[[1]]
  if (! identical(ai, afi)) stop("Different outcomes == molermat, molerfast")
  osize<-object.size(ac)
  tevs<-system.time(evs<-eigen(ac))[[1]]
  mtable[[ni,2]]<-ti
  mtable[[ni,3]]<-tc
  mtable[[ni,4]]<-osize
  mtable[[ni,5]]<-tevs
  mtable[[ni,6]]<-tfi
  mtable[[ni,7]]<-tfc
# cat(n, ti, tc, osize,"\n")
}

## ----molertimedisp, echo=FALSE-------------------------------------------
bmattym<-data.frame(n=mtable[,1], buildi=mtable[,2], buildc=mtable[,3], 
     osize=mtable[,4], eigentime=mtable[,5], bfast=mtable[,6],
     bfastc=mtable[,7])
print(bmattym)
cat("buildi - interpreted build time; buildc - byte compiled build time\n")
cat("osize - matrix size in bytes; eigentime - all eigensolutions time\n")
cat("bfast - interpreted vectorized build time; bfastc - same code, byte compiled time\n")

## ----drawtime1, echo=FALSE, fig.pos='h'----------------------------------
ti<-as.vector(mtable[,2])
tc<-as.vector(mtable[,3])
os<-as.vector(mtable[,4])
n<-as.vector(mtable[,1])
plot(n, ti)
title(main="Execution time vs matrix size")
title(sub="Regular Moler matrix routine, interpreted and byte compiled")
points(n, tc, pch=3, col='red')
legend(50,1,c("interpreted","byte compiled"), pch = c(1,3))
n2<-n*n
itime<-lm(ti~n+n2)
summary(itime)
ctime<-lm(tc~n+n2)
summary(ctime)
osize<-lm(os~n+n2)
summary(osize)

## ----rqdir, echo=TRUE----------------------------------------------------
rqdir<-function(x, AA){
  rq<-0.0
  n<-length(x) # assume x, AA conformable
  for (i in 1:n) {
     for (j in 1:n) {
        rq<-rq-x[i]*AA[[i,j]]*x[j] # Note - sign
     }
  }
  rq
}

## ----raynum1, echo=TRUE--------------------------------------------------
ray1<-function(x, AA){
    rq<- - t(x)%*%AA%*%x
}

## ----raynum2, echo=TRUE--------------------------------------------------
ray2<-function(x, AA){
    rq<- - as.numeric(crossprod(x, crossprod(AA,x)))
}

## ----raynum3, echo=TRUE--------------------------------------------------
ray3<-function(x, AA, ax=axftn){
    # ax is a function to form AA%*%x 
    rq<- - as.numeric(crossprod(x, ax(x, AA)))
}

## ----axm, echo=TRUE------------------------------------------------------
ax<-function(x, AA){
   u<- as.numeric(AA%*%x)
}

axx<-function(x, AA){
   u<- as.numeric(crossprod(AA, x))
}

## ----aximp, echo=TRUE----------------------------------------------------
aximp<-function(x, AA=1){ # implicit moler A*x
   n<-length(x)
   y<-rep(0,n)
   for (i in 1:n){
      tt<-0.
      for (j in 1:n) {
          if (i == j) tt<-tt+i*x[i]
          else tt<-tt+(min(i,j) - 2)*x[j]
      }
      y[i]<-tt 
   }
   y
}
ident<-function(x, B=1) x # identity

## ----axmfcode, echo=TRUE-------------------------------------------------
axmolerfast <- function(x, AA=1) {
# A fast and memory-saving version of A%*%x  
# For Moler matrix. Note we need a matrix argument to match other functions
n <- length(x)
j <- 1:n
ax <- rep(0, n)
for (i in 1:n) {
term <- x * (pmin(i, j) - 2)
ax[i] <- sum(term[-i]) 
}
ax <- ax + j*x
ax
}

## ----recompftn, echo=TRUE------------------------------------------------
system("rm moler.so")
system("rm moler.o")
system("R CMD SHLIB moler.f")
cat("Dynamic libraries rebuilt \n")

## ----axftn, echo=TRUE----------------------------------------------------
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")

axftn<-function(x, AA=1) { # ignore second argument
   n<-length(x) # could speed up by having this passed
   vout<-rep(0,n) # purely for storage
   res<-(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}

## ----cmpfns, echo=TRUE---------------------------------------------------
require(compiler)
axc<-cmpfun(ax)
axxc<-cmpfun(axx)
axftnc<-cmpfun(axftn)
aximpc<-cmpfun(aximp)
axmfc<-cmpfun(axmolerfast)

## ----timeax1, echo=TRUE--------------------------------------------------
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")
require(microbenchmark)
nmax<-10
ptable<-matrix(NA, nrow=nmax, ncol=11) # to hold results
# loop over sizes
for (ni in 1:nmax){
  n<-50*ni
  x<-runif(n) # generate a vector 
  ptable[[ni, 1]]<-n
  AA<-molermat(n)
  tax<-system.time(oax<-replicate(100,ax(x, AA))[,1])[[1]]
  taxc<-system.time(oaxc<-replicate(100,axc(x, AA))[,1])[[1]]
  if (! identical(oax, oaxc)) stop("oaxc NOT correct")
  taxx<-system.time(oaxx<-replicate(100,axx(x, AA))[,1])[[1]]
  if (! identical(oax, oaxx)) stop("oaxx NOT correct")
  taxxc<-system.time(oaxxc<-replicate(100,axxc(x, AA))[,1])[[1]]
  if (! identical(oax, oaxxc)) stop("oaxxc NOT correct")
  taxftn<-system.time(oaxftn<-replicate(100,axftn(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaxftn)) stop("oaxftn NOT correct")
  taxftnc<-system.time(oaxftnc<-replicate(100,axftnc(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaxftnc)) stop("oaxftnc NOT correct")
  taximp<-system.time(oaximp<-replicate(100,aximp(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaximp)) stop("oaximp NOT correct")
  taximpc<-system.time(oaximpc<-replicate(100,aximpc(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaximpc)) stop("oaximpc NOT correct")
  taxmfi<-system.time(oaxmfi<-replicate(100,axmolerfast(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaxmfi)) stop("oaxmfi NOT correct")
  taxmfc<-system.time(oaxmfc<-replicate(100,axmfc(x, AA=1))[,1])[[1]]
  if (! identical(oax, oaxmfc)) stop("oaxmfc NOT correct")
  ptable[[ni, 2]]<-tax
  ptable[[ni, 3]]<-taxc
  ptable[[ni, 4]]<-taxx
  ptable[[ni, 5]]<-taxxc
  ptable[[ni, 6]]<-taxftn
  ptable[[ni, 7]]<-taxftnc
  ptable[[ni, 8]]<-taximp
  ptable[[ni, 9]]<-taximpc
  ptable[[ni, 10]]<-taxmfi
  ptable[[ni, 11]]<-taxmfc
#  cat(n,tax, taxc, taxx, taxxc, taxftn, taxftnc, taximp, taximpc,"\n")
}
axtym<-data.frame(n=ptable[,1], ax=ptable[,2], axc=ptable[,3], 
  axx=ptable[,4], axxc=ptable[,5],
  axftn=ptable[,6], axftnc=ptable[,7], 
  aximp=ptable[,8], aximpc=ptable[,9], axmfast=ptable[,10],
  amfastc=ptable[,11])
print(axtym)

## ----adjaxtime, echo=TRUE------------------------------------------------
adjtym<-data.frame(n=axtym$n, axx1=axtym$axx+1*bmattym$buildi, 
     axxz=axtym$axx+100*bmattym$buildi, 
     axxc1=axtym$axxc+1*bmattym$buildc,axxcz=axtym$axxc+100*bmattym$buildc,
     axftn=axtym$axftn, aximp=axtym$aximp, aximpc=axtym$aximpc)
print(adjtym)

## ----rqtime1, echo=TRUE--------------------------------------------------
require(compiler)
rqdirc<-cmpfun(rqdir)
ray1c<-cmpfun(ray1)
ray2c<-cmpfun(ray2)
ray3c<-cmpfun(ray3)
dyn.load("moler.so")
  n<-500
  x<-runif(n) # generate a vector 
  AA<-molermat(n)
  tdi<-system.time(rdi<-replicate(100,rqdir(x, AA))[1])[[1]]
  tdc<-system.time(replicate(100,rdc<-rqdirc(x, AA))[1])[[1]]
  cat("Direct algorithm: interpreted=",tdi,"   byte-compiled=",tdc,"\n")
  t1i<-system.time(replicate(100,r1i<-ray1(x, AA))[1])[[1]]
  t1c<-system.time(replicate(100,r1c<-ray1c(x, AA))[1])[[1]]
  cat("ray1: mat-mult algorithm: interpreted=",t1i,"   byte-compiled=",t1c,"\n")
  t2i<-system.time(replicate(100,r2i<-ray2(x, AA))[1])[[1]]
  t2c<-system.time(replicate(100,r2c<-ray2c(x, AA))[1])[[1]]
  cat("ray2: crossprod algorithm: interpreted=",t2i,"   byte-compiled=",t2c,"\n")
  t3fi<-system.time(replicate(100,r3i<-ray3(x, AA, ax=axftn))[1])[[1]]
  t3fc<-system.time(replicate(100,r3i<-ray3c(x, AA, ax=axftnc))[1])[[1]]
  cat("ray3: ax Fortran + crossprod: interpreted=",t3fi,"   byte-compiled=",t3fc,"\n")
  t3ri<-system.time(replicate(100,r3i<-ray3(x, AA, ax=axmolerfast))[1])[[1]]
  t3rc<-system.time(replicate(100,r3i<-ray3c(x, AA, ax=axmfc))[1])[[1]]
  cat("ray3: ax fast R implicit + crossprod: interpreted=",t3ri,"   byte-compiled=",t3rc,"\n")

## ----rayspg1, echo=TRUE--------------------------------------------------
rqt<-function(x, AA){
    rq<-as.numeric(crossprod(x, crossprod(AA,x)))
}
proj<-function(x) { x/sqrt(crossprod(x)) }
require(BB)
n<-100
x<-rep(1,n)
AA<-molermat(n)
evs<-eigen(AA)
tmin<-system.time(amin<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=AA))[[1]]
#amin
tmax<-system.time(amax<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=-AA))[[1]]
#amax
evalmax<-evs$values[1]
evecmax<-evs$vectors[,1]
evecmax<-sign(evecmax[1])*evecmax/sqrt(as.numeric(crossprod(evecmax)))
emax<-list(evalmax=evalmax, evecmac=evecmax)
save(emax, file="temax.Rdata")
evalmin<-evs$values[n]
evecmin<-evs$vectors[,n]
evecmin<-sign(evecmin[1])*evecmin/sqrt(as.numeric(crossprod(evecmin)))
avecmax<-amax$par
avecmin<-amin$par
avecmax<-sign(avecmax[1])*avecmax/sqrt(as.numeric(crossprod(avecmax)))
avecmin<-sign(avecmin[1])*avecmin/sqrt(as.numeric(crossprod(avecmin)))
cat("minimal eigensolution: Value=",amin$value,"in time ",tmin,"\n")
cat("Eigenvalue - result from eigen=",amin$value-evalmin,"  vector max(abs(diff))=",
      max(abs(avecmin-evecmin)),"\n\n")
#print(amin$par)
cat("maximal eigensolution: Value=",-amax$value,"in time ",tmax,"\n")
cat("Eigenvalue - result from eigen=",-amax$value-evalmax,"  vector max(abs(diff))=",
      max(abs(avecmax-evecmax)),"\n\n")
#print(amax$par)

## ----runspg2, echo=FALSE-------------------------------------------------
require(compiler)
require(BB)
nmax<-10
stable<-matrix(NA, nrow=nmax, ncol=4) # to hold results
spgc<-cmpfun(spg)
rqtc<-cmpfun(rqt)
projc<-cmpfun(proj)

# loop over sizes
for (ni in 1:nmax){
  n<-50*ni
  x<-runif(n) # generate a vector 
  AA<-molerc(n) # make sure defined
  stable[[ni, 1]]<-n
  tbld<-system.time(AA<-molerc(n))[[1]]
  tspg<-system.time(aspg<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=-AA))[[1]]
  tspgc<-system.time(aspgc<-spgc(x, fn=rqtc, project=projc, control=list(trace=FALSE), AA=-AA))[[1]]
  stable[[ni, 2]]<-tspg
  stable[[ni, 3]]<-tspgc
  stable[[ni, 4]]<-tbld
#  cat(n,tspg, tspgc,tbld,"\n")
# times too short
}
spgtym<-data.frame(n=stable[,1], spgrqt=stable[,2], spgcrqtcaxc=stable[,3], tbldc=stable[,4])
print(spgtym)

## ----runopx1, echo=TRUE--------------------------------------------------
nobj<-function(x, AA=-AA){
   y<-x/sqrt(as.numeric(crossprod(x)))
   rq<- as.numeric(crossprod(y, crossprod(AA,y)))
}

ngrobj<-function(x, AA=-AA){
   y<-x/sqrt(as.numeric(crossprod(x))) 
   n<-length(x)
   dd<-sqrt(as.numeric(crossprod(x)))
   T1<-diag(rep(1,n))/dd
   T2<- x%o%x/(dd*dd*dd)
   gt<-T1-T2
   gy<- as.vector(2.*crossprod(AA,y))
   gg<-as.numeric(crossprod(gy, gt))
} 
require(optplus)
# mset<-c("L-BFGS-B", "BFGS", "CG", "spg", "ucminf", "nlm", "nlminb", "Rvmmin", "Rcgmin")
mset<-c("L-BFGS-B", "BFGS", "CG", "spg", "ucminf", "nlm", "nlminb", "Rvmmin", "Rcgmin")
nmax<-5
for (ni in 1:nmax){
  n<-20*ni
  x<-runif(n) # generate a vector 
  AA<-molerc(n) # make sure defined
  aall<-optimx(x, fn=nobj, gr=ngrobj, method=mset, AA=-AA, 
     control=list(starttests=FALSE, dowarn=FALSE))
  optansout(aall, NULL)
  cat("Above for n=",n," \n")
}

## ----rcgrun1, echo=TRUE--------------------------------------------------
ctable<-matrix(NA, nrow=10, ncol=2)
nmax<-10
for (ni in 1:nmax){
  n<-50*ni
  x<-runif(n) # generate a vector 
  AA<-molerc(n) # make sure defined
  tcgu<-system.time(arcgu<-Rcgminu(x, fn=nobj, gr=ngrobj, AA=-AA))[[1]]
  ctable[[ni,1]]<-n
  ctable[[ni,2]]<-tcgu
}
cgtime<-data.frame(n=ctable[,1], tRcgminu=ctable[,2])
print(cgtime)

## ----geradincode, echo=FALSE---------------------------------------------
ax<-function(x, AA){
   u<-as.numeric(AA%*%x)
}

bx<-function(x, BB){
   v<-as.numeric(BB%*%x)
}

geradin<-function(x, ax, bx, AA, BB, control=list(trace=TRUE, maxit=1000)){
# Geradin minimize Rayleigh Quotient, Nash CMN Alg 25
#  print(control)
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

## ----rungeradin1, echo=TRUE----------------------------------------------
cat("Test geradin with explicit matrix multiplication\n")
n<-10
AA<-molermat(n)
BB=diag(rep(1,n))
x<-runif(n)
tg<-system.time(ag<-geradin(x, ax, bx, AA=AA, BB=BB, 
   control=list(trace=FALSE)))[[1]]
cat("Minimal eigensolution\n")
print(ag)
cat("Geradin time=",tg,"\n")
tgn<-system.time(agn<-geradin(x, ax, bx, AA=-AA, BB=BB,
   control=list(trace=FALSE)))[[1]]
cat("Maximal eigensolution (negative matrix)\n")
print(agn)
cat("Geradin time=",tgn,"\n")

## ----timeger1, echo=TRUE-------------------------------------------------
naximp<-function(x, A=1){ # implicit moler A*x
   n<-length(x)
   y<-rep(0,n)
   for (i in 1:n){
      tt<-0.
      for (j in 1:n) {
          if (i == j) tt<-tt+i*x[i]
          else tt<-tt+(min(i,j) - 2)*x[j]
      }
      y[i]<- -tt # include negative sign
   }
   y
}

dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")

naxftn<-function(x, A) { # ignore second argument
   n<-length(x) # could speed up by having this passed
   vout<-rep(0,n) # purely for storage
   res<-(-1)*(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}

require(compiler)
naxftnc<-cmpfun(naxftn)
naximpc<-cmpfun(naximp)

require(microbenchmark)
nmax<-10
gtable<-matrix(NA, nrow=nmax, ncol=6) # to hold results
# loop over sizes
for (ni in 1:nmax){
  n<-50*ni
  x<-runif(n) # generate a vector 
  gtable[[ni, 1]]<-n
  AA<-molermat(n)
  BB<-diag(rep(1,n))
  tgax<-system.time(ogax<-geradin(x, ax, bx, AA=-AA, BB=BB, control=list(trace=FALSE)))[[1]]
  gtable[[ni, 2]]<-tgax
  tgaximp<-system.time(ogaximp<-geradin(x, naximp, ident, AA=1, BB=1, control=list(trace=FALSE)))[[1]]
  gtable[[ni, 3]]<-tgaximp
  tgaximpc<-system.time(ogaximpc<-geradin(x, naximpc, ident, AA=1, BB=1, control=list(trace=FALSE)))[[1]]
  gtable[[ni, 4]]<-tgaximpc
  tgaxftn<-system.time(ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)))[[1]]
  gtable[[ni, 5]]<-tgaxftn
  tgaxftnc<-system.time(ogaxftnc<-geradin(x, naxftnc, ident, AA=1, BB=1, control=list(trace=FALSE)))[[1]]
  gtable[[ni, 6]]<-tgaxftnc
#  cat(n,tgax, tgaximp, tgaximpc, tgaxftn, tgaxftnc,"\n")
}

gtym<-data.frame(n=gtable[,1], ax=gtable[,2], aximp=gtable[,3], 
  aximpc=gtable[,4], axftn=gtable[,5], axftnc=gtable[,6])
print(gtym)

## ----gercheck1, echo=TRUE------------------------------------------------
n<-100
x<-runif(n)
emax<-load("temax.Rdata")
evalmax<-emax$evalmax
evecmac<-emax$evecmax
ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE))
gvec<-ogaxftn$x
gval<- -ogaxftn$RQ
gvec<-sign(gvec[[1]])*gvec/sqrt(as.numeric(crossprod(gvec)))
diff<-gvec-evecmax
cat("Geradin diff eigenval from eigen result: ",gval-evalmax,"   max(abs(vector diff))=",
      max(abs(diff)), "\n")

## ----persp1, echo=FALSE, fig.pos='htbp'----------------------------------
cf<-data.frame(
    n=bmattym$n,
    eig=bmattym$eigentime,
    spg=spgtym$spgcrqtcaxc,
    rcg=cgtime$tRcgminu,
    ger=gtym$axftnc
)
eigen<-cf$eig/cf$ger
spg<-cf$spg/cf$ger
rcgmin<-cf$rcg/cf$ger
nsize<-cf$n
jn<-data.frame(nsize=nsize, eigen=eigen, spg=spg, rcgmin=rcgmin)
joe<-write.csv(jn, file="jndata.csv")
plot(nsize,spg, pch=1, xlab="n", ylab="time ratio")
points(nsize, rcgmin, pch=3)
points(nsize, eigen, pch=4)
title("Ratio of eigensolution times to Geradin routine by matrix size")
points(nsize, rep(1,10), type="l")
#legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4), lty = c(1,2,3))
legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4))


## ----n2000a, echo=FALSE--------------------------------------------------
dyn.load("moler.so")
n<-2000
t2000b<-system.time(AA<-molerc(n))[[1]]
t2000e<-system.time(evs<-eigen(AA))[[1]]
x<-runif(n)
t2000c<-system.time(ac<-Rcgminu(x, fn=nobj, gr=ngrobj, AA=-AA))[[1]]
t2000g<-system.time(ag<-geradin(x, naxftnc, ident, AA=1, BB=1, control=list(trace=FALSE)))[[1]]
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
cat("Ratios: build=", t2000b/t2000g, "eigen=",t2000e/t2000g,"  Rcgminu=",t2000c/t2000g,"\n")

## ----geradincodea, echo=TRUE---------------------------------------------
ax<-function(x, AA){
   u<-as.numeric(AA%*%x)
}
bx<-function(x, BB){
   v<-as.numeric(BB%*%x)
}
geradin<-function(x, ax, bx, AA, BB, control=list(trace=TRUE, maxit=1000)){
# Geradin minimize Rayleigh Quotient, Nash CMN Alg 25
# print(control)
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
  ans<-list(x=x, RQ=p0, ipr=ipr, msg=msg) # step 26
}

