f <- function(x) exp(x)*sin(3*x)*tanh(5*cos(30*x))
# fminbnd(f, -1, 1)
sfsp <- fminsfp(f, 0, 1)
print(sfsp)
