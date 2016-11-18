## @knitr runoptprob.R
runoptprob <- function(pfilename, minmeth=NULL, submeth=NULL, nstart=0, 
                 runopts=list(), control=list(), ...) {

  #- ?? Need to eval(parse()) ALL functions available, since f calls res etc.
  #- Need to carefully ensure these exist to avoid errors??
  #- ?? can we simplify and NOT have to eval(parse()) them, but simply source the prb file?
  
  print(runopts)
  print(control)
  optecho <- TRUE # temporarily at least, or put in a profile
  #- Get the path to the files (where should these be? Probably somehow related to pkg)   
  pfile <- paste(pfilename, ".prb", sep='')

  starts <- NA
  mformula <- NA # Make sure these are defined (they get set up in pfile)
#- ?? need to figure out dynamic setting of paths
  pfilepath <- paste("/home/john/rsvnall/optimizer/pkg/NISTO/inst/extdata/",pfile,sep='')
  source(pfilepath, echo=optecho) # -- filename (at least the root)
  cat("Objects in workspace:\n")
  print(ls())

  #- now have a lot of the information
  
  #-    -- which tool to use (optimr, nls, nlmrt, nlsr tools, minpack.LM tools)
  
  #-    -- choice of gradient function or approximation (gr= (gr, "grfwd", etc.))
  
  #-    -- controls -- as per the control list in programs
  
  #-    -- other arguments
  
  #-    -- xdata or dotargs (how to specify might be interesting)
  
  #-    -- timing control (e.g., microbenchmark or simple timing)
  
  #-  - read output control profile (initially just use sink())
  #-  -- make sure we have time/date stamp on all runs
  fname<-paste(pfilename, format(Sys.time(), "%Y%m%d%H%M"),".out",sep='')
  #- ?? not created until later, then conditionally 
  #- ?? sink(fname, append=TRUE, split=TRUE)
  
  #-  - read the file and execute it (make sure it has **R** commands so we can
  #-   actually source() it)
#  setupfn <- eval(parse(text=paste(pfilename,".setup", sep=''))) #- setup function
#  pdat <- setupfn()
#  cat("pdat:\n")
#  print(pdat)
  #- get the data and the starts
#  print(pdat$df)
#  dfname <- pdat$df
  starts <- eval(parse(text=paste(pfilename, ".starts", sep='')))
  uformula <-  eval(parse(text=paste(pfilename, ".formula", sep='')))
  udata <-  eval(parse(text=paste(pfilename, ".df", sep='')))
  ufn <- eval(parse(text=paste(pfilename,".f", sep='')))
  ures <- eval(parse(text=paste(pfilename,".res", sep='')))
  ujac <- eval(parse(text=paste(pfilename,".jac", sep='')))
  ugr <- eval(parse(text=paste(pfilename,".g", sep='')))

  havestarts <- exists(starts)
  haveuformula <- exists(uformula)
  cat("starts and formula: ", havestarts, haveuformula,"\n")

  #- - analyze the call to runprob and do the appropriate call
  
  #- - format output and extract and store summaries
  
  #-  -- this may be multilayerd and take a lot of work
  
  #-  -- start with no formatting, and gradually add features
  
  #-  -- need to save conditions
  if (nstart == 0) { # need to loop
    cat("nstart == 0 \n")
    print(starts)
    nst <- seq(1, (dim(starts)[1]))
  } else { nst <- nstart }
  if (minmeth == 'nls') {r
    for (istrt in nst){
      #- ?? need to extract options and arguments like trace
      #- some documentation output needed ??
      sol <- nls(uformula, data=udata, start=starts[istrt,], trace=TRUE)      
      print(sol)
      print(summary(sol))
    }    
  }  
  if (minmeth == 'nlxb') {
##??      require(nlmrt)
      for (istrt in nst){
      #- ?? need to extract options and arguments like trace
      #- some documentation output needed ??
      sol <- nlxb(uformula, data=udata, start=starts[istrt,], trace=TRUE)      
      print(sol)
      print(summary(sol))
    }    
  }  
  if (minmeth == "optimr") {
## ??     require(optimrx) #- ?? optimr for CRAN
    #- here need to check if they exist??
     if( is.null(runopts$gr) || ! is.character(runopts$gr) ) {
#       #- name.gr now a function      
       ugr <- eval(parse(text=paste(pfilename,".g", sep='')))
     } else { ugr <- (runopts$gr) }
     ufn <-  eval(parse(text=paste(pfilename,".f", sep='')))
     for (istrt in nst){
        sol <- optimr(starts[istrt,], ufn, ugr, method=submeth, control=list(trace=1))
        print(sol)
     }
  }
    #- result should be a list of things run
  #- ?? sink() # should make conditional and not usually do this, else delete after asking
  testsol <- list(pfilename=pfilename, nstart=nstart, minmeth=minmeth, submeth=submeth)
#- ?? definitely need more stuff returned
}