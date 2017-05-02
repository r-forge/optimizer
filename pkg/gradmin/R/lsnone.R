lsnone<-function(fn, fbest, xc, d, grv, ws, ...) { 
  # lsnone -- ALWAYS returns unit stop
  # ?? count fevals?
  rlout <- 1.0
  ws$nf <- ws$nf+1
  attr(rlout, "Fval")<- fn(xc + d)
  attr(rlout, "FAIL")<- FALSE # can never fail
  rlout
} # end lsnone
