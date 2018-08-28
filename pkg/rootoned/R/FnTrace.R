TraceSetup <- function(ifn=0, igr=0, ftrace=FALSE, fn=NA, gr=NA){
# JN: Define globals here
##    groot<-list(ifn=ifn, igr=igr, ftrace=ftrace, fn=fn, gr=gr, label="none")
##    envroot <<- list2env(groot) # Note globals in FnTrace
   ## if (getRversion() >= '2.15.1') utils::globalVariables(c('envroot'))
   ## utils::globalVariables("envroot") # Try declaring here -- causes errors
   ## This generates a NOTE that 
   ## TraceSetup: no visible binding for '<<-' assignment to ‘envroot’
   ## envroot<-list2env(groot, parent=.GlobalEnv) # Note globals in globals.R
# end globals
   envroot$ifn <- ifn
   envroot$igr <- igr
   envroot$ftrace <- ftrace
   envroot$fn <- fn
   envroot$gr <- gr
   return()
##   envroot
}

FnTrace <- function(x,...) { 
  # Substitute function to call when rootfinding
  # Evaluate fn(x, ...)
    val <- envroot$fn(x, ...)
    envroot$ifn <- envroot$ifn + 1 # probably more efficient ways
    if (envroot$ftrace) {
       cat("f(",x,")=",val," after ",envroot$ifn," ",envroot$label,"\n")
    }
    val
}

grTrace <- function(x,...) { 
  # Substitute function to call when rootfinding
  # Evaluate fn(x, ...)
  val <- envroot$gr(x, ...)
  envroot$igr <- envroot$igr + 1 # probably more efficient ways
  if (envroot$ftrace) {
    cat("gr(",x,")=",val," after ",envroot$igr," ",envroot$label,"\n")
  }
  val
}


FTTest <- function(){
  fsin <- function(x, fpar=0.5){ sin(x) - fpar }
  gsin <- function(x) { cos(x) }
  TraceSetup(ftrace=TRUE, fn=fsin, gr=gsin)
  val1 <- FnTrace(1)
  print(val1)
  ri <- c(0, 1.5)
  guess <- c(0,NA)
  #  print(str(envroot))
  tr <- rootwrap(fn=fsin, ri=ri, method="root1d")
  tr
  #  print(str(envroot))
  tu <- rootwrap(fn=fsin, ri=ri, method="uniroot")
  tu
#  print(str(envroot))
  tz <- rootwrap(fn=fsin, ri=ri, method="zeroin")
  tz
#  print(str(envroot))
  tb <- rootwrap(fn=fsin, ri=ri, method="bisect")
  tb
#  print(str(envroot))
  trf <- rootwrap(fn=fsin, ri=ri, method="regulaFalsi")
  trf
#  print(str(envroot))
  tm <- rootwrap(fn=fsin, ri=ri, method="muller")
  tm
#  print(str(envroot))
  tbr <- rootwrap(fn=fsin, ri=ri, method="brent")
  tbr
#  print(str(envroot))
  tn <- rootwrap(fn=fsin, ri=guess, method="newton")
  tn
#  print(str(envroot))
  tn1 <- rootwrap(fn=fsin, gr=gsin, ri=guess, method="newt1d")
  tn1
}  
