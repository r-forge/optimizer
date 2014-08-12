modeldoc <- function(fun, savefile=NULL) {
  efun <- environment(fun) # get the environment with data, expression, etc
  avn <- all.vars(efun$modelformula) # vars and parameters
  pnames <- names(efun$pvec) # assumes that vector is named (normal)
# ?? what do we do when it is not -- need to name it p1, p2, etc.
  iprm <- match(pnames, avn)
  funname <- match.call()[[2]]
  if (is.character(savefile)) { sink(savefile, split=TRUE) }
  else if (! is.null(savefile)) stop("savefile must be null or a character filename")
#  cat("Function name:",funname,"\n")
  notprm <- avn[-iprm] # get the non-parameter names
  dnames <- ls(efun$data) # names of vars in data dataframe
  idata <- match(dnames, notprm) # index for vars in data dataframe
  extradata <- notprm[-idata] # other data 
  cat("Content of function",funname,"\n")
  cat("The code:\n")
  print(modelexpr(fun))
  cat("Parameters:")
  print(pnames)
  cat("Data in 'data' environment:")
  print(dnames)
  cat("Extra data needed:")
  print(extradata)
  cat("Values:\n")
  cat("Parameters:")
  print(efun$pvec)
  cat("Data in dataframe:\n")
  for (vv in seq_along(dnames)){
     dident <- paste("efun$data$",dnames[vv],sep='')
     cat(dident,":")
     print(eval(parse(text=dident)))
  }
  cat("Extra data (!! NOT WORKING RIGHT !!):\n")
  if (length(extradata) < 1) cat("None\n")
  else {
    for (vv in seq_along(extradata)){
      cat(extradata[vv],":")
      print(eval(parse(text=vv)))
    }
  }
  if (is.character(savefile)) { sink() }
  return(0)
}
