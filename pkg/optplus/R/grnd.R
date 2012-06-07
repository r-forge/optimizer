############### grnd.R ####################
## require("numDeriv")
grnd<-function(par, userfn, ...) { # using grad from numDeriv
# with Richardson method
   tryg<-grad(userfn, par, ...)
}
############### end grnd ####################

