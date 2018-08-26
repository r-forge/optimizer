multrfind <- function(fn, gr=NULL, ri, ftrace=FALSE, meths="uniroot"){
## ?? no dot args??
  # run multiple methods for rootfinding
  TraceSetup(ftrace=ftrace)
  envroot$fn<-fn
  envroot$gr<-gr
  fguess <- c(ri[1],NA)
  #  allmeth <- c( "uniroot", "root1d", "zeroin", "newt1d", "newton", "bisect", "secant", "regulaFalsi", "muller", "brent")
  
  cnames <- c("root", "froot", "rtol", "iter", "fncount", "method")
  resdf <- as.data.frame(matrix(NA, nrow=length(meths), ncol=length(cnames)))
  colnames(resdf) <- cnames                         
  # print(str(resdf))
  irow <- 0
  for (method in meths) {
    irow <- irow + 1
    if (ftrace) {cat("Method ",irow," = ", method,"\n")}
    tri <- ri
    envroot$label <- method
    envroot$ifn <- 0
    envroot$igr <- 0
    # fix start for GUESS methods
    if (method %in% c("newt1d", "newton") ) { tri[2] <- NA }
    res <- rootwrap(fn, gr, ri, method=method, ftrace=ftrace)
    # ?? Store result
    # print(str(res))
    res <- res[cnames]
    resdf[irow,] <-  res
  } # end for loop
  ## Clean up
  resdf <- resdf[order(abs(resdf$froot)),]
  resdf
} 
