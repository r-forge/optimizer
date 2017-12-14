# To run these tests:
#   library(optimx)
#   library(svUnit)
#   runit.all <- system.file("unit1", "runit.1.R", package = "optimx")
#   source(runit.1); clearLog(); test.1()
#   Log()

test.1 <- function() {

	## test 1

	checkTrue(require(optimx))

	XRosenbrock.f <- function (x) 
	{
	    n <- length(x)
	    x1 <- x[2:n]
	    x2 <- x[1:(n - 1)]
	    sum(100 * (x1 - x2^2)^2 + (1 - x2)^2)
	}
	XRosenbrock.g <- function (x) 
	{
	    n <- length(x)
	    g <- rep(NA, n)
	    g[1] <- 2 * (x[1] - 1) + 400 * x[1] * (x[1]^2 - x[2])
	    if (n > 2) {
	        ii <- 2:(n - 1)
	        g[ii] <- 2 * (x[ii] - 1) + 400 * x[ii] * (x[ii]^2 - x[ii + 
	            1]) + 200 * (x[ii] - x[ii - 1]^2)
	    }
	    g[n] <- 200 * (x[n] - x[n - 1]^2)
	    g
	}

	xstart <- c(-1.2,1)

	ansone<-optimr(xstart, XRosenbrock.f, XRosenbrock.g, method="Rvmmin")
#	ansone.sum<-summary(ansone, order= "value")[1,]
  tt <- str(ansone)



        ## test 3
	ansoneb<-optimx(c(1), f1, lower=c(-1), upper=c(10),control=list(all.methods=TRUE))
	ansoneb.sum<-summary(ansoneb, order= "value")[1,]
	ansoneb.sum.target<- structure(list(p1 = 0.555555555557684, value = 2.35390946502058, 
	    fevals = 17, gevals = NA_real_, niter = NA_real_, convcode = 0, 
	    kkt1 = TRUE, kkt2 = TRUE, xtimes = 0), .Names = c("p1", "value", 
	    "fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes"
	    ), details = structure(list("bobyqa", 5.52926311436595e-11, 
	    structure(14.0000000000441, .Dim = c(1L, 1L)), 14.0000000000441, "none"), 
	    .Dim = c(1L, 5L), .Dimnames = list("bobyqa", c("method", "ngatend", "nhatend",
	    "hev", "message" ))), maximize = FALSE, npar = 1L, row.names = "bobyqa", 
	    class = c("optimx", "data.frame"))


	# don't compare xtimes
	ansoneb.sum$xtimes <- ansoneb.sum.target$xtimes <- NULL
	checkEquals(ansoneb.sum, ansoneb.sum.target)


}


ansone.target <- structure(list(par = c(1, 1), value = 0, counts = structure(c(59, 
39), .Names = c("function", "gradient")), convergence = 2, message = "Rvmminu appears to have converged"), .Names = c("par", 
"value", "counts", "convergence", "message"))
print(checkEquals(ansone, ansone.target))


