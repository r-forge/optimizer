wrapnls <-function(formula, start, trace=FALSE, data=NULL, lower=NULL, upper=NULL, masked=NULL, control=list(),...){
#
#  A wrapper to call nlsmnq() and then call nls() with the solution.
#  The calling sequence matches that of nlsmnq()
#
    if (is.null(data)) stop("wrapnls() must have 'data' supplied")
    if (!is.null(masked)) stop("nls() cannot handle masked parameters")
    if (is.null(lower) && is.null(upper)) bounded<-FALSE else bounded<-TRUE
    pnames<-names(start)
    if (bounded) { # deal with bounds constrained nonlinear least squares
       first<-nlsmnqb(formula, start, trace=trace, data=data, lower=lower,
                       upper=upper, control=control,...)
       #  Should check this has worked, but ...
       newstart<-first$coeffs
       names(newstart)<-pnames
       print(first)
       #  Should put this in a try(), but let's see it work first
       second<-nls(formula, newstart, trace=trace, data=data, 
              algorithm='port', lower=lower, upper=upper, control=control,...)
     } else {
       first<-nlsmnq(formula, start, trace=trace, data=data, control=control,...)
       #  Should check this has worked, but ...
       newstart<-first$coeffs
       names(newstart)<-pnames
       print(first)
       #  Should put this in a try(), but let's see it work first
       second<-nls(formula, newstart, trace=trace, data=data, control=control,...)
     }
}
