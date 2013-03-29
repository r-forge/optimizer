
# To run these tests:
#   library(optimx)
#   library(svUnit)
#   runit.all <- system.file("unitTests", "runit.all.R", package = "optimx")
#   source(runit.all); clearLog(); test.all()
#   Log()

test.all <- function() {

	## test 1

	checkTrue(require(optimx))

	## test 2

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
		return(gg)
	}

	genrose.h <- function(x, gs=NULL) { ## compute Hessian
	   if(is.null(gs)) { gs=100.0 }
		n <- length(x)
		hh<-matrix(rep(0, n*n),n,n)
		for (i in 2:n) {
			z1<-x[i]-x[i-1]*x[i-1]
			z2<-1.0-x[i]
			hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
			hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
			hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
			hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
		}
		return(hh)
	}

	startx<-4*seq(1:10)/3.
	ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
	print(ans8)
	print(ans8[, "gevals"])
	print(ans8["spg", ])
	print(ans8, par.select = 1:3)
	print(ans8, best.only = TRUE)

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
		return(gg)
	}

	genrose.h <- function(x, gs=NULL) { ## compute Hessian
	   if(is.null(gs)) { gs=100.0 }
		n <- length(x)
		hh<-matrix(rep(0, n*n),n,n)
		for (i in 2:n) {
			z1<-x[i]-x[i-1]*x[i-1]
			z2<-1.0-x[i]
			hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
			hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
			hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
			hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
		}
		return(hh)
	}

	startx<-4*seq(1:10)/3.
	ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
	ans8.sum <- summary(ans8, order = value)[1, ]
	ans8.sum.target <- structure(list(p1 = 0.999999999920414, 
	p2 = 0.999999999913547, 
	p3 = 0.999999999885138, p4 = 0.999999999895512, 
        p5 = 0.999999999967326, 
	p6 = 1.00000000007184, p7 = 1.00000000015134, 
        p8 = 1.00000000020547, 
	p9 = 1.00000000022618, p10 = 1.00000000021113, value = 1, 
	fevals = 147, gevals = 85, niter = NA_real_, convcode = 0, 
	kkt1 = 1, kkt2 = NA_real_, xtimes = 0.0199999999999996), 
        .Names = c("p1", 
	"p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "value", 
	"fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes"
	), details = structure(list(c(-4.15795776335918, -5.96304183914467, 
	23.5408226631421, 138.69710944406, -135.573412666897, 
        -66.6938432145199, 
	145.749092585896, -149.848543104669, 430.230152331936, 36.0982006928723
	), structure(c(-26.2301378857652, -5.94101708907067, 0, 0, 0, 
	0, 0, 0, 0, 0, -5.94101708907067, 36.816607304302, -28.8773141897125, 
	0, 0, 0, 0, 0, 0, 0, 0, -28.8773141897125, 144.030417517769, 
	-47.7258383065507, 0, 0, 0, 0, 0, 0, 0, 0, -47.7258383065507, 
	257.737719650414, -48.8012556369575, 0, 0, 0, 0, 0, 0, 0, 0, 
	-48.8012556369575, 235.828858468625, 57.1205282696378, 0, 0, 
	0, 0, 0, 0, 0, 0, 57.1205282696378, 16.6987694987266, -30.8777477665626, 
	0, 0, 0, 0, 0, 0, 0, 0, -30.8777477665626, 378.130491934927, 
	-76.8088785364331, 0, 0, 0, 0, 0, 0, 0, 0, -76.8088785364331, 
	338.396911215821, -86.3397947169126, 0, 0, 0, 0, 0, 0, 0, 0, 
	-86.3397947169126, 3031.64183470906, -242.695100166075, 0, 0, 
	0, 0, 0, 0, 0, 0, -242.695100166075, 22), .Dim = c(10L, 10L)), 
	c(3053.81693062519, 438.091835147444, 309.786714826021, 
        278.648277449347, 218.315287956348, 131.725821133951, 
        29.4580546463029, 2.41158302879628, 
	    -0.369729705835123, -26.8333026936891), "none"), .Dim = c(1L, 
	4L), .Dimnames = list("Nelder-Mead", c("ngatend", "nhatend", 
	"hev", "message"))), maximize = FALSE, npar = 10L, 
        row.names = "Rvmmin", class = c("optimx", "data.frame"))

	# don't compare xtimes
	ans8.sum$xtimes <- ans8.sum.target$xtimes <- NULL
	checkEquals(ans8.sum, ans8.sum.target)

}

