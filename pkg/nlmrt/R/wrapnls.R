wrapnls <-function(formula, start, trace=FALSE, data=NULL, control=list(),...){
#
#  A wrapper to call nlsmnq() and then call nls() with the solution.
#  The calling sequence matches that of nlsmnq()
#
    if (is.null(data)) stop("wrapnls() must have 'data' supplied")

    first<-nlsmnq(formula, start, trace=trace, data=data, control=control,...)
#  Should check this has worked, but ...
    newstart<-first$coeffs
    print(first)
#  Should put this in a try(), but let's see it work first
    second<-nls(formula, newstart, trace=trace, data=data, control=control,...)
}
