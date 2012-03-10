rm(list=ls()) # This one close, but no cigar.
   Form2resfun <- function(f, p = quote(p)) {
        xx <- all.vars(f)
        fp <- match(p, xx) # Problem in matching the names of params
        xx <- c(xx[fp], xx[-fp])
        ff <- vector("list", length(xx))
        names(ff) <- xx
        sf<-as.character(f)
        lhs<-sf[2] # NOTE ORDER formula with ~ puts ~, lhs, rhs
        rhs<-sf[3]
# And build the residual at the parameters
        resexp<-paste(rhs,"-",lhs, collapse=" ")
        ff[[length(ff) + 1]] <- crossprod(eval(parse(text=resexp)))
#  want crossprod(resexp)
        as.function(ff, parent.frame())
    }
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
    t<-1:length(y) # for testing
    f<- y ~ b[1]/(1+b[2]*exp(-b[3]*t))
    p<-c(b1=1, b2=1, b3=1)
    b<-p
    john<-Form2resfun(f, p)
    john(t)
