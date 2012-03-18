rm(list=ls()) # This one close, but no cigar.
#   Form2resfun <- function(f, p = quote(p)) {
    Form2resfun <- function(f, p ) {
        cat("In Form2resfun\n")
        xx <- all.vars(f)
        cat("f:")
        print(f)
        cat("p:")
        print(p)
        cat("xx:")
        print(xx)
        fp <- match(names(p), xx) # Problem in matching the names of params
        cat("fp:")
        print(fp)
        xx2 <- c(xx[fp], xx[-fp])
        ff <- vector("list", length(xx2))
        names(ff) <- xx2
        cat("ff:")
        print(ff)
        sf<-as.character(f)
        if ((length(sf)!=3) && (sf[1]!="~")) stop("Bad model formula expression")
        lhs<-sf[2] # NOTE ORDER formula with ~ puts ~, lhs, rhs
        rhs<-sf[3]
# And build the residual at the parameters
        resexp<-paste(rhs,"-",lhs, collapse=" ")
        cat("resexp:")
        print(resexp)
        fnexp<-paste("crossprod(",resexp,")", sep="")
        ff[[length(ff) + 1]] <- parse(text=fnexp)
#  want crossprod(resexp)
        myfn<-as.function(ff, parent.frame())
        cat("final fn:")
        print(myfn)
        myfn
    }
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
    t<-1:length(y) # for testing
#    f<- y ~ b[1]/(1+b[2]*exp(-1*b[3]*t))
    f<- y ~ b1/(1+b2*exp(-1*b3*t))
    p<-c(b1=1, b2=1, b3=1)
    b<-p
    john<-Form2resfun(f, p)
    ans<-eval(john(t=t,y=y, p))
    ans

