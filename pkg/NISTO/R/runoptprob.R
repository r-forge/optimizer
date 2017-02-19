## @knitr runoptprob.R
runoptprob <- function(pfilename, probclass=NULL, minmeth=NULL, submeth=NULL, istart="1", 
                 runopts=list(), control=list(), ...) {

  # Have to figure out WHAT we want to do
  #   - for a single problem -- apply all possible methods
  #   - this takes multiple passes / logic. May want to simplify.
  # How to proceed?
  #    What do we do about bounds?
  #      - always impose if present
  #      - have 2 problems if can omit bounds
  # Need probclass to decide between using squared .res vs. .f form of problem,
  #   but these may be the same in some cases. Could/should leave out .f, .g in such cases,
  #   and let this program sort it out.
  # 
  #- now have a lot of (the) information
    #-    -- which tool to use (optimr, nls, nlmrt, nlsr tools, minpack.LM tools)
    #-    -- choice of gradient function or approximation (gr= (gr, "grfwd", etc.))
    #-    -- controls -- as per the control list in programs
    #-    -- other arguments
    #-    -- xdata or dotargs (how to specify might be interesting)
    #-    -- timing control (e.g., microbenchmark or simple timing)
  #- ?? Need to eval(parse()) ALL functions available, since f calls res etc.
  #- Need to carefully ensure these exist to avoid errors??
  #- ?? can we simplify and NOT have to eval(parse()) them, but simply source the prb file?

  #- ?? Do we need to rm() all the things we test for, namely,
  #-   lower, upper, ...
  #- ?? should make bounds (probname).lower etc.
  
  solveformula <- c("stats::nls", "nlmrt::nlxb", "minpack.lm::nls.lm", "nls2::nls2")

  solvesumsquares <- c("nlmrt::nlfb", "minpack.lm::nlsLM")

  if (minmeth == "optimr") minmeth <- "optimrx" # update to optimrx to get more submeths

  solveuncopt <- c("optimrx::Nelder-Mead",
                   "optimrx::BFGS", 
                   "optimrx::CG",
                   "optimrx::hjn",
                   "optimrx::ucminf",
                   "optimrx::lbfgs",
                   "optimrx::spg",
                   "optimrx::uobyqa",
                   "optimrx::nlm",
                   "optimrx::subplex")

  solveboundopt <- c("optimrx::L-BFGS-B",
                   "optimrx::Rvmmin",
                   "optimrx::Rcgmin",
                   "optimrx::Rtnmin",
                   "optimrx::hjn",
                   "optimrx::hjkb",
                   "optimrx::nlminb",
                   "optimrx::nmkb",
                   "optimrx::lbfgsb3",
                   "optimrx::bobyqa",
                   "optimrx::nlminb")
#- Diagnostics  
  print(runopts)
  print(control)
  optecho <- TRUE # temporarily at least, or put in a profile
  #- Get the path to the files (where should these be? Probably somehow related to pkg)   
  #- Set up counts
  if (! exists("counters")) { counters <- new.env() }
  counters$kf <- 0
  counters$kg <- 0
  counters$kres <- 0
  counters$kjac <- 0
  counters$khess <- 0
  counters$kform <- 0 # How to use this ??
  #- end counts setup
  
  pfile <- paste(pfilename, ".prb", sep='')

  starts <- NA
  mformula <- NA # Make sure these are defined (they get set up in pfile)
#- ?? need to figure out dynamic setting of paths
  ## ?? how to change to a well-defined path based on where pkg installed??
  pfilepath <- paste("/home/john/rsvnall/optimizer/pkg/NISTO/inst/extdata/",pfile,sep='')
  source(pfilepath, echo=optecho) # -- filename (at least the root)
  cat("Objects in workspace:\n")
  print(ls())

  #-  - read output control profile (initially just use sink())
  #-  -- make sure we have time/date stamp on all runs
  fname<-paste(pfilename, format(Sys.time(), "%Y%m%d%H%M"),".out",sep='')
  #- ?? not created until later, then conditionally 
  #- ?? sink(fname, append=TRUE, split=TRUE)
  
  #-  - read the file and execute it (make sure it has **R** commands so we can
  #-   actually source() it)
  cstarts <- paste(pfilename, ".start", sep='')
  cuformula <-  paste(pfilename, ".formula", sep='')
  cudata <-  paste(pfilename, ".df", sep='')
  cufn <- paste(pfilename,".f", sep='')
  cures <- paste(pfilename,".res", sep='')
  cujac <- paste(pfilename,".jac", sep='')
  cugr <- paste(pfilename,".g", sep='')
  culower <- paste(pfilename,".lower", sep='')
  cuupper <- paste(pfilename,".upper", sep='')
  havestarts <- FALSE
  if (exists(cstarts)) {
     havestarts <- TRUE
     print(cstarts)
     fstart <- eval(parse(text=cstarts))
     print(fstart)
     print(istart)
     strt <- fstart(istart) #- Push names into the probfile
  }
  haveuformula <- FALSE
  if (exists(cuformula)) {
    uformula <-  eval(parse(text=cuformula))
    haveuformula <- TRUE
  }
  haveudata <- FALSE
  if (exists(cudata)) {
    udata <-  eval(parse(text=cudata))
    haveudata <- TRUE
  }
  haveufn <- FALSE
  if (exists(cufn)) {
    ufn <- eval(parse(text=cufn))
    haveufn <- TRUE
  }
  haveugr <- FALSE
  if (exists(cugr)) {
    ugr <- eval(parse(text=cugr))
    haveugr <- TRUE
  }
  haveures <- FALSE
  if (exists(cures)) {
    ures <- eval(parse(text=cures))
    haveures <- TRUE
  }
  haveujac <- FALSE
  if (exists(cujac)) {
    ujac <- eval(parse(text=cujac))
    haveujac <- TRUE
  }

  havebounds <- (exists(culower) || exists(cuupper))
  #- Do not need to parse -- these are already either parsed or don't exist
  
  #- classes of problems
  pclass = c("uncopt", "sumsquares", "formula", "boundopt")  
  upclass <- c() #- ?? may not need this
  #- Work out all possible tools for this PROBLEM file (disregard what call requests)
  #- Could put this in a separate file of available tools in form method::submeth
  if (havebounds) {
    if (haveuformula) {
      solveform <- c("stats::nls", "nlmrt::nlxb", "minpack.lm::nls.lm")
      #- nls2 does not mention bounds
    }
    if (haveures){
      #- ?? build function from res, gradient from jac. 
      #- ?? need specified and/or default derivative approach?
      solveopt <- solveboundopt
    }  
    if (haveufn){
      solveopt <- solveboundopt
      #- use function, set gradient from code and/or approximations
    }
  } else {
      if (haveuformula) {
          solveform <- solveformula
      } 
      upclass <- c(upclass, "formula")
      if (haveures){
        #- ?? build function from res, gradient from jac. 
        #- ?? need specified and/or default derivative approach?
        solveopt <- solveuncopt
      }  
      if (haveufn){
        solveopt <- solveuncopt
        #- use function, set gradient from code and/or approximations
      }
  } # end set up list of solvers
    
    upclass <- c(upclass, "boundopt") 
  if (haveufn){ upclass <- c(upclass, "uncopt")}
  if (haveures) {upclass <- c(upclass, "sumsquares")}
  
  #- - analyze the call to runprob and do the appropriate call
  #- - format output and extract and store summaries
  #-  -- this may be multilayerd and take a lot of work
  #-  -- start with no formatting, and gradually add features
  #-  -- need to save conditions

  if (minmeth == 'nls') {
      sol <- nls(uformula, data=udata, start=strt, trace=TRUE)      
      print(sol)
      print(summary(sol))
  }  
  if (minmeth == 'nlxb') {
##??      require(nlmrt)
      #- ?? need to extract options and arguments like trace
      #- some documentation output needed ??
      sol <- nlxb(uformula, data=udata, start=strt, trace=TRUE)      
      print(sol)
      print(summary(sol))
  }  
  if (minmeth == "optimrx") {
## ??     require(optimrx) #- ?? optimr for CRAN
    #- here need to check if they exist??
     if( is.null(runopts$gr) || ! is.character(runopts$gr) ) {
#       #- name.gr now a function      
       ugr <- eval(parse(text=paste(pfilename,".g", sep='')))
     } else { ugr <- (runopts$gr) }
     ufn <-  eval(parse(text=paste(pfilename,".f", sep='')))
     sol <- optimr(strt, ufn, ugr, method=submeth, control=list(trace=1))
     print(sol)
   }
    #- result should be a list of things run
  #- ?? sink() # should make conditional and not usually do this, else delete after asking
  testsol <- list(pfilename=pfilename, istart=istart, minmeth=minmeth, submeth=submeth)
#- ?? definitely need more stuff returned
}
