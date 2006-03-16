if(!require("numDeriv"))stop("this test requires numDeriv.")
Sys.info()
  

#######################################################################

# Test gradient and hessian calculation in genD using  data for calculating 
#   curvatures in Bates and Watts.
#model A p329,data set 3 (table A1.3, p269) Bates & Watts (Puromycin example)

#######################################################################
  


puromycin <- function(th){
   x <- c(0.02,0.02,0.06,0.06,0.11,0.11,0.22,0.22,0.56,0.56,1.10,1.10)
   y <- c(76,47,97,107,123,139,159,152,191,201,207,200)
   ( (th[1] * x)/(th[2] + x) ) - y
  }

D.anal <- function(th){
  # analytic derivatives. Note numerical approximation gives a very good
  # estimate of these, but neither give D below exactly. The results are very
  # sensitive to th, so rounding error in the reported value of th could explain
  # the difference. But more likely th is correct and D has been rounded for
  # publication - and the analytic D with published th seems to work best.
  # th = c(212.70188549 ,  0.06410027) is the nls est of th for BW published D.
  x <- c(0.02,0.02,0.06,0.06,0.11,0.11,0.22,0.22,0.56,0.56,1.10,1.10)
  y <- c(76,47,97,107,123,139,159,152,191,201,207,200)
  cbind(x/(th[2]+x), -th[1]*x/(th[2]+x)^2, 0, -x/(th[2]+x)^2, 2*th[1]*x/(th[2]+x)^3)
 }

# D matrix from p235. This may be useful for rough comparisons, but rounding
# used for publication introduces substantial errors. check D.anal1 - D.BW
D.BW <- t(matrix(c(
0.237812, -601.458, 0, -2.82773, 14303.4,
0.237812, -601.458, 0, -2.82773, 14303.4,
0.483481, -828.658, 0, -3.89590, 13354.7,
0.483481, -828.658, 0, -3.89590, 13354.7,
0.631821, -771.903, 0, -3.62907, 8867.4,
0.631821, -771.903, 0, -3.62907, 8867.4,
0.774375, -579.759, 0, -2.72571, 4081.4,
0.774375, -579.759, 0, -2.72571, 4081.4,
0.897292, -305.807, 0, -1.43774,  980.0,
0.897292, -305.807, 0, -1.43774,  980.0,
0.944936, -172.655, 0, -0.81173,  296.6,
0.944936, -172.655, 0, -0.81173,  296.6),   5,12))

D.anal <- D.anal(c(212.7000, 0.0641))

D.calc <- genD(puromycin,c(212.7000, 0.0641)) # compares to 1e-4 below
# increasing r does not always help
#D.calc <- genD(puromycin,c(212.7000, 0.0641), r=10)#compares to 0.01 below
#D.calc <- genD(puromycin,c(212.7000, 0.0641), d=0.001)

cat("numerical D:")
print( D.calc$D, digits=16)

cat("diff. between analytic and numerical D:")
print( D.calc$D - D.anal, digits=16)

cat("max. abs. diff. between analtic and numerical D:")
print( max(abs(D.calc$D-D.anal)), digits=16)

# These would be interesting except for 0 column
#cat("% diff. between analtic and numerical D:")
#z <- 100 * (D.calc$D - D.anal) / D.anal
#print( z, digits=16)

#cat("max. abs. % diff. between analtic and numerical D:")
#print( max(abs(z)), digits=16)

if(max(abs(D.calc$D - D.anal)) < 1e-4) invisible(T) else stop("BW test FAILED")

