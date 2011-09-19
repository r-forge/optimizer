#see page 857 or page xxx (30) of vol II of Facchinei & Pang (2003)

#Fischer-Burmeister
phiFB <- function(a, b) 
	sqrt(a^2+b^2) - (a+b)
GrAphiFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, a / sqrt(a^2+b^2) - 1)
GrBphiFB <- function(a, b) 
	ifelse(a == 0 & b == 0, -1/2, b / sqrt(a^2+b^2) - 1)

#minimum
phiMin <- function(a, b) 
	min(a, b)
GrAphiMin <- function(a, b) 
	1*(a <= b)
GrBphiMin <- function(a, b) 
	1*(b <= a)

#Mangasarian
phiMan <- function(a, b, f)
	f(abs(a-b)) - f(a) - f(b)

GrAphiMan <- function(a, b, fprime)
	sign(a-b) * fprime(abs(a-b)) - fprime(a)

GrBphiMan <- function(a, b, fprime)
	sign(b-a) * fprime(abs(a-b)) - fprime(b)

#Luo-Tseng
phiLT <- function(a, b, q)
	(a^q + b^q)^(1/q) - (a+b)

GrAphiLT <- function(a, b, q)
	ifelse(a == 0 & b == 0,
		1 + (1/2)^((q-1)/q),
		1 - sign(a)*( abs(a) / (a^q + b^q)^(1/q) )^(q-1) )

GrBphiLT <- function(a, b, q)
	ifelse(a == 0 & b == 0,
		1 + (1/2)^((q-1)/q),
		1 - sign(b)*( abs(b) / (a^q + b^q)^(1/q) )^(q-1) )



#Kanzow-Kleinmichel
phiKK <- function(a, b, lambda)
	(sqrt( (a-b)^2 + 2*lambda*a*b ) - (a+b) ) / (2-lambda)

GrAphiKK <- function(a, b, lambda)
	ifelse( a == 0 & b == 0, 
		-1+sqrt(16*(2-lambda)*(lambda^2 - 2*lambda + 4))/8/(2-lambda),
		( (a+b*(lambda-1)) / sqrt( (a-b)^2 + 2*lambda*a*b ) - 1 ) / (2-lambda) )

GrBphiKK <- function(a, b, lambda)
	ifelse( a == 0 & b == 0, 
		   -1+sqrt(16*(2-lambda)*(lambda^2 - 2*lambda + 4))/8/(2-lambda),
		   ( (b+a*(lambda-1)) / sqrt( (a-b)^2 + 2*lambda*a*b ) - 1 ) / (2-lambda) )



