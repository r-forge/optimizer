
#Fischer-Burmeister
phiFB <- function(a, b) 
	sqrt(a^2+b^2) - (a+b)
GrAphiFB <- function(a, b) 
	a / sqrt(a^2+b^2) - 1
GrBphiFB <- function(a, b) 
	b / sqrt(a^2+b^2) - 1

#minimum
phiMin <- function(a, b) 
	min(a, b)
GrAphiMin <- function(a, b) 
	1*(a <= b)
GrBphiMin <- function(a, b) 
	1*(b <= a)

#Mangasarian
phiMan <- function(a, b, f=function(t) t^3 )
	f(abs(a-b)) - f(a) - f(b)

#Luo-Tseng
phiLT <- function(a, b, q=2)
	(sum( (a-b)^q ))^(1/q) - (a+b)

#page xxx or 857 of vol II