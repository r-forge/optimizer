# A one-dimensional ‘wild‘ function, taken from an example in the DEoptim package. 
# It has many local minima, and a global minimum of approx. -15.81515. 
#
# Karl Ove Hufthammer 2009/05/06 07:28

wild.f <- function(x) 10 * sin(0.3 * x) * sin(1.3 * x^2) + 1e-05 * x^4 + 0.2 * x + 80
wild.g <- function(x) 3 * cos(0.3 * x) * sin(1.3 * x^2) + 26 * sin(0.3 * x) * cos(1.3 * 
  x^2) * x + 4e-05 * x^3 + 0.2
wild.h <- function(x) -0.9 * sin(0.3 * x) * sin(1.3 * x^2) + 15.6 * cos(0.3 * x) * 
  cos(1.3 * x^2) * x - 67.6 * sin(0.3 * x) * sin(1.3 * x^2) * 
  x^2 + 26 * sin(0.3 * x) * cos(1.3 * x^2) + 0.00012 * x^2

 
curve(wild.f, -50, 50, n = 1000)
