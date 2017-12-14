require(optimx)

#####################
## All examples are in this .Rd file
##
## Rosenbrock Banana function
fr <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
ansrosenbrock <- optimr(fn=fr,gr="grfwd", par=c(1,2), method="Rvmmin")
print(ansrosenbrock) # use print to allow copy to separate file that 
cat("No gr specified as a test\n")
ansrosenbrock0 <- optimr(fn=fr, par=c(1,2), method="Rvmmin")
print(ansrosenbrock0) # use print to allow copy to separate file that 
#    can be called using source()
#####################
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}
genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	gg
}

# analytic gradient test
xx<-rep(pi,10)
lower<-NULL
upper<-NULL
bdmsk<-NULL
genrosea<-optimr(xx,genrose.f, genrose.g, method="Rvmmin", gs=10)
genrosenf<-optimr(xx,genrose.f, gr="grfwd", method="Rvmmin", gs=10) # use local numerical gradient
genrosenullgr<-optimr(xx,genrose.f, method="Rvmmin", gs=10) # no gradient specified
cat("genrosea uses analytic gradient\n")
print(genrosea)
cat("genrosenf uses grfwd standard numerical gradient\n")
print(genrosenf)
cat("genrosenullgr has no gradient specified\n")
print(genrosenullgr)
cat("If optextras is loaded, then other numerical gradients can be used.\n")

cat("timings direct call Bounded vs Unbounded\n")
require(Rvmmin)
xx<-rep(pi,50)
lo<-rep(-100,50)
up<-rep(100,50)
bdmsk<-rep(1,50)
tb<-system.time(ab<-Rvmminb(xx,genrose.f, genrose.g, lower=lo, upper=up, bdmsk=bdmsk))[1]
tu<-system.time(au<-Rvmminu(xx,genrose.f, genrose.g))[1]
cat("times U=",tu,"   B=",tb,"\n")
cat("solution Rvmminu\n")
print(au)
cat("solution Rvmminb\n")
print(ab)
cat("diff fu-fb=",au$value-ab$value,"\n")
cat("max abs parameter diff = ", max(abs(au$par-ab$par)),"\n")
# Note: there are multiple solutions!

cat("timings optimr call Bounded vs Unbounded\n")
require(Rvmmin)
lo<-rep(-100,50)
up<-rep(100,50)
tuo<-system.time(auo<-optimr(xx,genrose.f, genrose.g, method="Rvmmin"))[1]
tbo<-system.time(abo<-optimr(xx,genrose.f, genrose.g, lower=lo, upper=up, method="Rvmmin"))[1]
cat("times U=",tuo,"   B=",tbo,"\n")
cat("solution optimr:Rvmmin unconstrained\n")
print(auo)
cat("solution optimr:Rvmmin bounded\n")
print(abo)
cat("diff fu-fb=",auo$value-abo$value,"\n")
cat("max abs parameter diff = ", max(abs(auo$par-abo$par)),"\n")

#####################
cat("test bounds and masks\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds1<-optimr(startx,genrose.f, genrose.g, lower=lo,upper=up, method="Rvmmin") 
print(grbds1)

cat("test lower bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
grbds2<-optimr(startx,genrose.f, genrose.g, lower=lo, method="Rvmmin") 
print(grbds2)

cat("test lower bound single value only\n")
nn<-4
startx<-rep(pi,nn)
lo<-2
up<-rep(10,nn)
grbds3<-optimr(startx,genrose.f, genrose.g, lower=lo, method="Rvmmin") 
print(grbds3)

cat("test upper bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds4<-optimr(startx,genrose.f, genrose.g, upper=up, method="Rvmmin") 
print(grbds4)

cat("test upper bound single value only\n")
nn<-4
startx<-rep(pi,nn)
grbds5<-optimr(startx,genrose.f, genrose.g, upper=10, method="Rvmmin") 
print(grbds5)



cat("test masks only\n")
nn<-6
bd<-c(1,1,0,0,1,1)
startx<-rep(pi,nn)
grbds6<-optimr(startx,genrose.f, genrose.g, bdmsk=bd, method="Rvmmin") 
print(grbds6)

cat("test upper bound on first two elements only\n")
nn<-4
startx<-rep(pi,nn)
upper<-c(10,8, Inf, Inf)
grbds7<-optimr(startx,genrose.f, genrose.g, upper=upper, method="Rvmmin") 
print(grbds7)


cat("test lower bound on first two elements only\n")
nn<-4
startx<-rep(0,nn)
lower<-c(0,1.1, -Inf, -Inf)
grbds8<-optimr(startx,genrose.f,genrose.g,lower=lower, method="Rvmmin", control=list(maxit=2000)) 
print(grbds8)

cat("test n=1 problem using simple squares of parameter\n")

sqtst<-function(xx) {
   res<-sum((xx-2)*(xx-2))
}

nn<-1
startx<-rep(0,nn)
onepar<-optimr(startx,sqtst, gr="grfwd", method="Rvmmin", control=list(trace=1)) 
print(onepar)

cat("Suppress warnings\n")
oneparnw<-optimr(startx,sqtst, gr="grfwd", method="Rvmmin", control=list(dowarn=FALSE,trace=1)) 
print(oneparnw)

#####################
grall <- opm(xx, genrose.f, genrose.g, method="ALL")
summary(grall, order=value)

