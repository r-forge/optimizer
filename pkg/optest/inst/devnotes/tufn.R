         tufn <- function(par, ...){
            f <- fr(par, ...)
            g <- grr(par, ...)
            attr(f,"gradient") <- g
            attr(f,"hessian") <- NULL # ?? maybe change later
            f
        }