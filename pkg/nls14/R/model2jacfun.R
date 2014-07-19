model2jacfun <- function(modelformula, pvec, funname = "myjac", 
    filename = NULL) {
    pnames <- names(pvec)
    # Creates Jacobian function
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
    ff <- vector("list", length(xx2))
    names(ff) <- xx2
    parts <- strsplit(as.character(es), "~")[[1]]
    if (length(parts) != 2) 
        stop("Model expression is incorrect!")
    lhs <- parts[1]
    rhs <- parts[2]
    # And build the residual at the parameters
    if (lhs == "") { # we allow 1-sided expressions 140716, but drop ~ for jacobian
       resexp <- paste(rhs, collapse = " ")
    } else {
       resexp <- paste(rhs, "-", lhs, collapse = " ")
    }
    jacexp <- deriv(parse(text = resexp), pnames)  # gradient expression
    jfstr <- paste("jstruc<-with(data,eval(", jacexp,"))", sep = "")  ##3
#?    cat("jfstr:")
#?    print(jfstr)
    pparse <- ""
    for (i in 1:npar) {
        pparse <- paste(pparse, "   ", pnames[[i]], "<-prm[[", 
            i, "]]\n", sep = "")
    }
     myjstr <- paste(funname, "<-function(prm, data=NULL) {\n", 
# Want NULL to trip error if we forget to call properly with data set.
        pparse, jfstr, " \n", "jacmat<-attr(jstruc,'gradient')\n ", 
        "return(jacmat)\n }", sep = "")
    if (!is.null(filename)) 
        write(myjstr, file = filename)  # write out the file
    myparse <- parse(text = myjstr)
    # This may be cause trouble if there are errors
    if (class(myparse) == "try-error") 
        stop("Error in Jacobian code string")
    eval(myparse)
}
