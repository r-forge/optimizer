model2ssfun <- function(modelformula, pvec, funname = "myss", 
    filename = NULL) {
    pnames <- names(pvec)
    if (is.null(pnames)) 
        stop("MUST have named parameters in pvec")
    if (is.character(modelformula)) {
        es <- modelformula
    } else {
        tstr <- as.character(modelformula)  # note ordering of terms!
        if (length(tstr) == 2) { # 1-sided formula
             es <- paste("~", tstr[[2]], "")
	} else {
	     es <- paste(tstr[[2]], "~", tstr[[3]], sep="")
        }
    }
    xx <- all.vars(parse(text = es))
    rp <- match(pnames, xx)  # Match names to parameters
    xx2 <- c(xx[rp], xx[-rp])
    xxparm <- xx[rp]
    npar <- length(xxparm)
    xxvars <- xx[-rp]
    nvar <- length(xxvars)
    parts <- strsplit(as.character(es), "~")[[1]]
    if (length(parts) != 2) 
        stop("Model expression is incorrect!")
    lhs <- parts[1]
    rhs <- parts[2]
    # And build the residual at the parameters
    if (lhs == "") { # we allow 1-sided expressions 140716
       resexp <- paste(rhs, collapse = " ")
    } else {
       resexp <- paste(rhs, "-", lhs, collapse = " ")
    }
    fnexp <- paste("resids<-as.numeric(eval(", resexp, "))", 
        sep = "")  ##3
    pparse <- ""
    for (i in 1:npar) {
        pparse <- paste(pparse, "   ", pnames[[i]], "<-prm[[", 
            i, "]]\n", sep = "")
    }
    myfstr <- paste(funname, "<-function(prm, data=NULL) {\n", 
        pparse, fnexp, "\n  ss<-as.numeric(crossprod(resids))\n }\n", 
        sep = "")
    if (!is.null(filename)) 
        write(myfstr, file = filename)  # write out the file
    tparse <- try(parse(text = myfstr))
    # This may cause trouble if there are errors
    if (class(tparse) == "try-error") 
        stop("Error in residual code string")
    eval(tparse)
}
