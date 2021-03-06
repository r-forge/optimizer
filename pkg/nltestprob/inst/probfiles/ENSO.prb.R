## @knitr ##ENSO.prb
# This is file ##ENSO.prb

rm(list=ls()) # May want to do this on each problem for safety.

probname <- "##ENSO"
probdesc <- "
NIST/ITL StRD
Dataset Name:  ENSO              (ENSO.dat)

File Format:   ASCII
Starting Values   (lines 41 to  49)
Certified Values  (lines 41 to  54)
Data              (lines 61 to 228)

Procedure:     Nonlinear Least Squares Regression

Description:   The data are monthly averaged atmospheric pressure 
differences between Easter Island and Darwin, 
Australia.  This difference drives the trade winds in 
the southern hemisphere.  Fourier analysis of the data
reveals 3 significant cycles.  The annual cycle is the
strongest, but cycles with periods of approximately 44
and 26 months are also present.  These cycles
correspond to the El Nino and the Southern Oscillation.
Arguments to the SIN and COS functions are in radians.

Reference:     Kahaner, D., C. Moler, and S. Nash, (1989). 
Numerical Methods and Software.  
Englewood Cliffs, NJ: Prentice Hall, pp. 441-445.

Data:          1 Response  (y = atmospheric pressure)
1 Predictor (x = time)
168 Observations
Average Level of Difficulty
Observed Data

Model:         Miscellaneous Class
9 Parameters (b1 to b9)

y = b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
+ b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
+ b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 )  + e

Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   11.0        10.0          1.0510749193E+01  1.7488832467E-01
b2 =    3.0         3.0          3.0762128085E+00  2.4310052139E-01
b3 =    0.5         0.5          5.3280138227E-01  2.4354686618E-01
b4 =   40.0        44.0          4.4311088700E+01  9.4408025976E-01
b5 =   -0.7        -1.5         -1.6231428586E+00  2.8078369611E-01
b6 =   -1.3         0.5          5.2554493756E-01  4.8073701119E-01
b7 =   25.0        26.0          2.6887614440E+01  4.1612939130E-01
b8 =   -0.3        -0.1          2.1232288488E-01  5.1460022911E-01
b9 =    1.4         1.5          1.4966870418E+00  2.5434468893E-01

Residual Sum of Squares:                    7.8853978668E+02
Residual Standard Deviation:                2.2269642403E+00
Degrees of Freedom:                               159
Number of Observations:                           168





Data:   y          x
12.90000    1.000000
11.30000    2.000000
10.60000    3.000000
11.20000    4.000000
10.90000    5.000000
7.500000    6.000000
7.700000    7.000000
11.70000    8.000000
12.90000    9.000000
14.30000   10.000000
10.90000    11.00000
13.70000    12.00000
17.10000    13.00000
14.00000    14.00000
15.30000    15.00000
8.500000    16.00000
5.700000    17.00000
5.500000    18.00000
7.600000    19.00000
8.600000    20.00000
7.300000    21.00000
7.600000    22.00000
12.70000    23.00000
11.00000    24.00000
12.70000    25.00000
12.90000    26.00000
13.00000    27.00000
10.90000    28.00000
10.400000    29.00000
10.200000    30.00000
8.000000    31.00000
10.90000    32.00000
13.60000    33.00000
10.500000    34.00000
9.200000    35.00000
12.40000    36.00000
12.70000    37.00000
13.30000    38.00000
10.100000    39.00000
7.800000    40.00000
4.800000    41.00000
3.000000    42.00000
2.500000    43.00000
6.300000    44.00000
9.700000    45.00000
11.60000    46.00000
8.600000    47.00000
12.40000    48.00000
10.500000    49.00000
13.30000    50.00000
10.400000    51.00000
8.100000    52.00000
3.700000    53.00000
10.70000    54.00000
5.100000    55.00000
10.400000    56.00000
10.90000    57.00000
11.70000    58.00000
11.40000    59.00000
13.70000    60.00000
14.10000    61.00000
14.00000    62.00000
12.50000    63.00000
6.300000    64.00000
9.600000    65.00000
11.70000    66.00000
5.000000    67.00000
10.80000    68.00000
12.70000    69.00000
10.80000    70.00000
11.80000    71.00000
12.60000    72.00000
15.70000    73.00000
12.60000    74.00000
14.80000    75.00000
7.800000    76.00000
7.100000    77.00000
11.20000    78.00000
8.100000    79.00000
6.400000    80.00000
5.200000    81.00000
12.00000    82.00000
10.200000    83.00000
12.70000    84.00000
10.200000    85.00000
14.70000    86.00000
12.20000    87.00000
7.100000    88.00000
5.700000    89.00000
6.700000    90.00000
3.900000    91.00000
8.500000    92.00000
8.300000    93.00000
10.80000    94.00000
16.70000    95.00000
12.60000    96.00000
12.50000    97.00000
12.50000    98.00000
9.800000    99.00000
7.200000   100.00000
4.100000   101.00000
10.60000   102.00000
10.100000   103.00000
10.100000   104.00000
11.90000   105.00000
13.60000    106.0000
16.30000    107.0000
17.60000    108.0000
15.50000    109.0000
16.00000    110.0000
15.20000    111.0000
11.20000    112.0000
14.30000    113.0000
14.50000    114.0000
8.500000    115.0000
12.00000    116.0000
12.70000    117.0000
11.30000    118.0000
14.50000    119.0000
15.10000    120.0000
10.400000    121.0000
11.50000    122.0000
13.40000    123.0000
7.500000    124.0000
0.6000000    125.0000
0.3000000    126.0000
5.500000    127.0000
5.000000    128.0000
4.600000    129.0000
8.200000    130.0000
9.900000    131.0000
9.200000    132.0000
12.50000    133.0000
10.90000    134.0000
9.900000    135.0000
8.900000    136.0000
7.600000    137.0000
9.500000    138.0000
8.400000    139.0000
10.70000    140.0000
13.60000    141.0000
13.70000    142.0000
13.70000    143.0000
16.50000    144.0000
16.80000    145.0000
17.10000    146.0000
15.40000    147.0000
9.500000    148.0000
6.100000    149.0000
10.100000    150.0000
9.300000    151.0000
5.300000    152.0000
11.20000    153.0000
16.60000    154.0000
15.60000    155.0000
12.00000    156.0000
11.50000    157.0000
8.600000    158.0000
13.80000    159.0000
8.700000    160.0000
8.600000    161.0000
8.600000    162.0000
8.700000    163.0000
12.80000    164.0000
13.20000    165.0000
14.00000    166.0000
13.40000    167.0000
14.80000    168.0000

"

#- Note: environment / list "counters" must already exist

if (exists("pe")) { 
  rm("pe")  
}

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

#- nls format expression
ENSO.formula <- ( y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) +
                              b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 ) +
                              b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ) )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("ENSO")))# Optimization test function enso
# enso from NISTnls
# ??ref...


enso.f <- function(x) {
   res<-enso.res(x)
   f<-sum(res*res)
}

enso.res <- function(b) {
# NOTE: could benefit from some sort of constraint to avoid equal parameters in trig args.
   xx<-ENSO$x # note case!
   yy<-ENSO$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   b9<-b[9]
   res<-b1 + b2*cos( 2*pi*xx/12 ) + b3*sin( 2*pi*xx/12 ) + b5*cos( 2*pi*xx/b4 ) + b6*sin( 2*pi*xx/b4 ) + b8*cos( 2*pi*xx/b7 ) + b9*sin( 2*pi*xx/b7 )  - yy
   return(res)
}

# enso - Jacobian
enso.jac <- function(b) {
stop("not defined")
   xx<-enso$x
   yy<-enso$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

enso.h <- function(x) {
stop("not defined")
   JJ<-enso.jac(x)
   H <- t(JJ) %*% JJ
   res<-enso.res(x)
stop("not defined")

}

enso.g<-function(x) {
#   stop("not defined")
   JJ<-enso.jac(x)
   res<-enso.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

enso.fgh<-function(x) {
   f<-enso.f(x)
   g<-enso.g(x)
   H<-enso.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1 <- c( 11.0,  3.0, 0.5, 40.0, -0.7, -1.3, 25.0, -0.3, 1.4)
names(start1) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9")
start2 <- c(10.0, 3.0, 0.5, 44.0, -1.5, 0.5, 26.0, -0.1, 1.5) 
names(start2) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9")
       
## Examples

cat("Examples from NISTnls:\n")

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = ENSO)
plot(y ~ x, data = ENSO, type = "l")  # to see the pattern more clearly
Try(fm1 <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
               + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
               + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
               data = ENSO, trace = TRUE,
               start = c(b1 = 11.0, b2 = 3.0, b3 = 0.5, b4 = 40.0, b5 = -0.7,
                         b6 = -1.3, b7 = 25.0, b8 = -0.3, b9 = 1.4)))
Try(fm1a <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
                data = ENSO, trace = TRUE, alg = "port",
                start = c(b1 = 11.0, b2 = 3.0, b3 = 0.5, b4 = 40.0, b5 = -0.7,
                          b6 = -1.3, b7 = 25.0, b8 = -0.3, b9 = 1.4)))
Try(fm2 <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
               + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
               + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
               data = ENSO, trace = TRUE,
               start = c(b1 = 10.0, b2 =  3.0, b3 =  0.5, b4 = 44.0, b5 = -1.5,
                         b6 =  0.5, b7 = 26.0, b8 = -0.1, b9 =  1.5)))
Try(fm2a <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
                data = ENSO, trace = TRUE, alg = "port",
                start = c(b1 = 10.0, b2 =  3.0, b3 =  0.5, b4 = 44.0, b5 = -1.5,
                          b6 =  0.5, b7 = 26.0, b8 = -0.1, b9 =  1.5)))
Try(fm3 <- nls(y ~ cbind(1, cos( 2*pi*x/12 ), sin( 2*pi*x/12 ), cos( 2*pi*x/b4 ),
                         sin( 2*pi*x/b4 ), cos( 2*pi*x/b7 ), sin( 2*pi*x/b7 )),
               data = ENSO, trace = TRUE,
               start = c(b4 = 40.0, b7 = 25.0), algorithm = "plinear"))
Try(fm4 <- nls(y ~ cbind(1, cos( 2*pi*x/12 ), sin( 2*pi*x/12 ), cos( 2*pi*x/b4 ),
                         sin( 2*pi*x/b4 ), cos( 2*pi*x/b7 ), sin( 2*pi*x/b7 )),
               data = ENSO, trace = TRUE,
               start = c(b4 = 44.0, b7 = 26.0), algorithm = "plinear"))


library(nlsr)
ENSO1nlxb <- nlxb(start=start1, trace=TRUE, formula=ENSO.formula, data=mypdata)
print(ENSO1nlxb)

ENSO2nlxb <- nlxb(start=start2, trace=TRUE, formula=ENSO.formula, data=mypdata)
print(ENSO2nlxb)

library(optimr)
ENSO1opmfwd <- opm(start1, enso.f, "grfwd", method="ALL")
summary(ENSO1opmfwd, order=value)
ENSO2opmfwd <- opm(start2, enso.f, "grfwd", method="ALL")
summary(ENSO2opmfwd, order=value)
