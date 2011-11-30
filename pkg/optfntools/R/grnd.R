############### grnd.R ####################
## require("numDeriv")
grnd<-function(par, userfn, ...) { # using grad from numDeriv
#grnd<-function(par, ...) { # using grad from numDeriv
##   userfn<-attr(grnd,"userfn") #?? test
##   cat("grnd userfn: ")
##   print(userfn)
   tryg<-grad(userfn, par, ...)
}
############### end grnd ####################

