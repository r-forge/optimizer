ssums<-function(x){
  n<-length(x)
  tt<-sum(x)
  ss<-1:n
  xx<-(x/tt)*(x/tt)
  sum(ss*xx)
}

cat("Try penalized sum\n")
require(optimx)
st<-runif(3)
aos<-optimx(st, ssums, control=list(all.methods=TRUE))
# rescale the parameters
nsol<-dim(aos)[1]
for (i in 1:nsol){ 
#  tpar<-aos[[i,2]] # second part of each solution is the parameter vector
  tpar<-aos[[i,1]] # first column of each solution is the parameter vector
  aos[[i,1]]<-tpar/sum(tpar)
}
optansout(aos, NULL)

# Implementing linear equality constraint in `spg'
# sum(x) = 1
#  Note: a direct projection approach: x/sum(x) does not work
#

ssum<-function(x){
  n<-length(x)
  ss<-1:n
  xx<-x*x
  sum(ss*xx)
}

Amat <- matrix(c(1,1,1), 1, 3, byrow=TRUE)
b <- 1
meq <- 1  # one equality constraint:  Ax = b

spg(par=st, fn=ssum, project="projectLinear",  projectArgs=list(A=Amat, b=b, meq=meq))

