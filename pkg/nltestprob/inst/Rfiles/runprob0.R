runoptprob <- function(pfilename=NULL, minmeth='nls', submeth=NULL, 
                       nstart = 0, options=list(args=NULL, control=NULL) ) {
  #- ?? Need to eval(parse()) ALL functions available, since f calls res etc.
  #- Need to carefully ensure these exist to avoid errors??
  #- ?? can we simplify and NOT have to eval(parse()) them, but simply source the prb file?
  
  print(options$args)
  print(options$controls)
  optecho <- TRUE # temporarily at least, or put in a profile
  #- Get the path to the files (where should these be? Probably somehow related to pkg)   
  pfile <- paste(pfilename, ".prb", sep='')
  source(pfile, echo=optecho) # -- filename (at least the root)
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
  sink(fname, append=TRUE, split=TRUE)
  
  #-  - read the file and execute it (make sure it has **R** commands so we can
  #-   actually source() it)
  eval(parse(text=psetup.txt)) #- setup function
  #- get the data and the starts
  pds <- eval(parse(text=paste(probname,".setup()",sep='')))
  
         
  #- - analyze the call to runprob and do the appropriate call
  
  #- - format output and extract and store summaries
  
  #-  -- this may be multilayerd and take a lot of work
  
  #-  -- start with no formatting, and gradually add features
  
  #-  -- need to save conditions
  if (nstart == 0) { # need to loop
    nst <- seq(1: dim(pds$starts)[1])
  } else { nst <- nstart }
  if (minmeth == 'nls') {
    for (istrt in nst){
      #- ?? need to extract options and arguments like trace
      #- some documentation output needed ??
      sol <- nls(mformula, data=pds$df, start=pds$starts[istrt,], trace=TRUE)      
      print(sol)
      print(summary(sol))
    }    
  }  
  if (minmeth == 'nlxb') {
      require(nlmrt)
      for (istrt in nst){
      #- ?? need to extract options and arguments like trace
      #- some documentation output needed ??
      sol <- nlxb(mformula, data=pds$df, start=pds$starts[istrt,], trace=TRUE)      
      print(sol)
      print(summary(sol))
    }    
  }  
  if (minmeth == "optimr") {
     require(optimrx)
     eval(parse(text=fn.txt))
     if( is.null(options$args[["gr"]]) || ! is.character(options$args[["gr"]]) ) {
       eval(parse(text=gr.txt))
       #- name.gr now a function      
       ugr <- eval(parse(text=paste(probname,".g", sep='')))
     } else { ugr <- (options$args[["gr"]]) }
     ufn <-  eval(parse(text=paste(probname,".f", sep='')))
     for (istrt in nst){
        sol <- optimr(pds$starts[istrt,], ufn, ugr, method=submeth, control=list(trace=1))
        print(sol)
     }
  }
    #- result should be a list of things run
  sink()
  testsol <- list(probname=probname, nstart=nstart, minmeth=minmeth, submeth=submeth)
}

# This attempts to call up a test problem and run it

pfname <- "DanielWood"
# try it
test1 <- runoptprob(pfilename=pfname, minmeth="nls")
test1
tmp <- readline("continue")

test2 <- runoptprob(pfilename=pfname, minmeth="nlxb")
test2
tmp <- readline("continue")

test3 <- runoptprob(pfilename=pfname, minmeth="optimr", submeth="Rvmmin")
test3

