## @knitr runoptprob.R
runoptprob <- function(pfilename, minmeth=NULL, submeth=NULL, nstart=0, 
                 runopts=list(), control=list(), ...) {

  #- ?? Need to eval(parse()) ALL functions available, since f calls res etc.
  #- Need to carefully ensure these exist to avoid errors??
  #- ?? can we simplify and NOT have to eval(parse()) them, but simply source the prb file?

  solveformula <- c("stats::nls", "nlmrt::nlxb", "minpack.lm::nls.lm")
  solvesumsquares <- c("nlmrt::nlfb", "minpack.lm::nlsLM")
  if (minmeth == "optimr") minmeth <- "optimrx" # update to optimrx to get more submeths
  solveuncopt <- c("optimrx::Nelder-Mead",
                   "optimrx::BFGS", 
                   "optimrx::L-BFGS-B",
                   "optimrx::CG",
                   "optimrx::Rvmmin",
                   "optimrx::Rcgmin",
                   "optimrx::Rtnmin",
                   "optimrx::hjn",
                   "optimrx::hjkb",
                   "optimrx::nmkb",
                   "optimrx::ucminf",
                   "optimrx::lbfgsb3",
                   "optimrx::lbfgs",
                   "optimrx::spg",
                   "optimrx::bobyqa",
                   "optimrx::uobyqa",
                   "optimrx::nlm",
                   "optimrx::nlminb")

  solveboundopt <- c("optimrx::L-BFGS-B",
                   "optimrx::Rvmmin",
                   "optimrx::Rcgmin",
                   "optimrx::Rtnmin",
                   "optimrx::hjn",
                   "optimrx::hjkb",
                   "optimrx::nmkb",
                   "optimrx::lbfgsb3",
                   "optimrx::bobyqa",
                   "optimrx::nlminb")
  
  print(runopts)
  print(control)
  optecho <- TRUE # temporarily at least, or put in a profile
  #- Get the path to the files (where should these be? Probably somehow related to pkg)   
  pfile <- paste(pfilename, ".prb", sep='')

  starts <- NA
  mformula <- NA # Make sure these are defined (they get set up in pfile)
#- ?? need to figure out dynamic setting of paths
  ## ?? how to change to a well-defined path based on where pkg installed??
  pfilepath <- paste("/home/john/rsvnall/optimizer/pkg/NISTO/inst/extdata/",pfile,sep='')
  source(pfilepath, echo=optecho) # -- filename (at least the root)
  cat("Objects in workspace:\n")
  print(ls())

# Have to figure out WHAT we want to do
#   - for a single problem -- apply all possible methods
#   - this takes multiple passes / logic. May want to simplify.
# How to proceed?
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

  #- now have a lot of (the) information
  
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
  cstarts <- paste(pfilename, ".starts", sep='')
  cuformula <-  paste(pfilename, ".formula", sep='')
  cudata <-  paste(pfilename, ".df", sep='')
  cufn <- paste(pfilename,".f", sep='')
  cures <- paste(pfilename,".res", sep='')
  cujac <- paste(pfilename,".jac", sep='')
  cugr <- paste(pfilename,".g", sep='')
  havestarts <- FALSE
  if (exists(cstarts)) {
     starts <- eval(parse(text=cstarts))
     havestarts <- TRUE
  }
  
  haveuformula <- FALSE
  if (exists(cstarts)) {
    uformula <-  eval(parse(text=cuformula))
    haveuformula <- TRUE
  }
  
  haveudata <- FALSE
  if (exists(cstarts)) {
    udata <-  eval(parse(text=cudata))
    haveudata <- TRUE
  }
  
  haveufn <- FALSE
  if (exists(cstarts)) {
    ufn <- eval(parse(text=cufn))
    haveufn <- TRUE
  }
  
  haveures <- FALSE
  if (exists(cstarts)) {
    ures <- eval(parse(text=cures))
    haveures <- TRUE
  }
  
  haveujac <- FALSE
  if (exists(cstarts)) {
    ujac <- eval(parse(text=cujac))
    haveujac <- TRUE
  }
  
  haveugr <- FALSE
  if (exists(cstarts)) {
    ugr <- eval(parse(text=cugr))
    haveugr <- TRUE
  }

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
  if (minmeth == 'nls') {
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
  if (minmeth == "optimrx") {
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
