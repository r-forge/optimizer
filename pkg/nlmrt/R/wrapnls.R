wrapnls <-function(formula, start, trace=FALSE, data, lower=NULL, upper=NULL,
        control=list(),...){
#
#  A wrapper to call nlsmnq() and then call nls() with the solution.
#  The calling sequence matches that of nlsmnq()
#
    if (is.null(data)) stop("wrapnls() must have 'data' supplied")
# Note that there are no bounds or masks.
    first<-nlxb(formula, start, trace=trace, data=data, lower=lower, upper=upper, control=control,...)
#  Should check this has worked, but ...
    newstart<-first$coeffs
    print(first)
#  Should put this in a try(), but let's see it work first
    if ((is.null(lower) && is.null(upper)) || (all(lower)==(-Inf)) && (all(upper) == Inf)) {
       second<-nls(formula, newstart, trace=trace, data=data, control=control,...)
    } else {
       second<-nls(formula, newstart, trace=trace, data=data, algorithm="port",
            lower=lower, upper=upper, control=control,...)
    }
}
