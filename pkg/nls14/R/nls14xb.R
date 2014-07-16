nls14xb <- function(formula, start, trace = FALSE, data = NULL, 
    lower = -Inf, upper = Inf, masked = NULL, 
    weights=NULL, control = list(), 
    ...) {
    # A simplified and hopefully robust alternative to finding
    # the nonlinear least squares minimizer that causes
    # 'formula' to give a minimal residual sum of squares.
    # 
    # nls14xb is particularly intended to allow for the
    # resolution of very ill-conditioned or else near
    # zero-residual problems for which the regular nls()
    # function is ill-suited. 
    # 
    # J C Nash 2014-7-16   nashjc _at_ uottawa.ca
    # 
    # formula looks like 'y~b1/(1+b2*exp(-b3*T))' start MUST be
    # a vector where all the elements are named: e.g.,
    # start=c(b1=200, b2=50, b3=0.3) trace -- TRUE for console
    # output data is a data frame containing data for variables
    # used in the formula that are NOT the parameters. This
    # data may also be defined in the parent frame i.e.,
    # 'global' to this function lower is a vector of lower
    # bounds upper is a vector of upper bounds masked is a
    # character vector of names of parameters that are fixed.
    # control is a list of control parameters. These are: ...
    # 
    # ... will need to contain data for other variables that
    # appear in the formula and are defined in a parent frame
    # (Not sure how needed??) ?? need to fix.
    # 
    # This variant uses a qr solution without forming the sum
    # of squares and cross products t(J)%*%J
    # 
    # Function to display SS and point
    showpoint <- function(SS, pnum) {
        pnames <- names(pnum)
        npar <- length(pnum)
        cat("lamda:", lamda, " SS=", SS, " at")
        for (i in 1:npar) {
            cat(" ", pnames[i], "=", pnum[i])
        }
        cat(" ", feval, "/", jeval)
        cat("\n")
    }
# ?? need to sort out and maybe build a dataframe.
# ?? get names of any data in args or ...
# ?? No data, then create frame
# ?? else if data in args, add to data frame
# ?? or should we make user do this?
# ?? and put in the weights
    ######### get data from data frame if exists
    ######### print(str(data))
    if (!is.null(data)) {
        for (dfn in names(data)) {
            cmd <- paste(dfn, "<-data$", dfn, "")
            eval(parse(text = cmd))
        }
    } else stop("'data' must be a list or an environment")
    # ensure params in vector
    pnames <- names(start)
    start <- as.numeric(start)
    names(start) <- pnames
    # bounds
    npar <- length(start)  # number of parameters
    if (length(lower) == 1) 
        lower <- rep(lower, npar)
    if (length(upper) == 1) 
        upper <- rep(upper, npar)
    # ?? more tests on bounds
    if (length(lower) != npar) 
        stop("Wrong length: lower")
    if (length(upper) != npar) 
        stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) 
        stop("Infeasible start")
    if (trace) {
        cat("formula: ")
        print(formula)
        cat("lower:")
        print(lower)
        cat("upper:")
        print(upper)
    }
    # controls
    ctrl <- list(watch = FALSE, phi = 1, lamda = 1e-04, offset = 100, 
        laminc = 10, lamdec = 4, femax = 10000, jemax = 5000, rofftest = TRUE, 
        smallsstest = TRUE)
     ##   maxlamda <- 1e+60) ## dropped 130709 ??why?
    epstol <- (.Machine$double.eps) * ctrl$offset
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
        if (!(onename %in% nctrl)) {
            if (trace) 
                cat("control ", onename, " is not in default set\n")
        }
        ctrl[onename] <- control[onename]
    }
    if (trace) 
        print(ctrl)
    phiroot <- sqrt(ctrl$phi)
    lamda <- ctrl$lamda # Note spelling -- a throwback to Ag Can 1974 and
	## way to see if folk are copying code.
    offset <- ctrl$offset
    laminc <- ctrl$laminc
    lamdec <- ctrl$lamdec  # save typing
    watch <- ctrl$watch
    femax <- ctrl$femax
    jemax <- ctrl$jemax
    # First get all the variable names:
    vn <- all.vars(parse(text = formula))
    # Then see which ones are parameters (get their positions
    # in the set xx
    pnum <- start  # may simplify later??
    pnames <- names(pnum)
    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(pnames %in% masked)  # NOTE: %in% not == or order gives trouble
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters
    if (trace) {
        parpos <- match(pnames, vn)
        datvar <- vn[-parpos]  # NOT the parameters
        for (dvn in datvar) {
            cat("Data variable ", dvn, ":")
            print(eval(parse(text = dvn)))
        }
    }
    if (is.character(formula)) {
        es <- formula
    } else {
        tstr <- as.character(formula)  # note ordering of terms!
        es <- paste(tstr[[2]], "~", tstr[[3]], "")
    }
    # Now separate the sides
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
    for (i in 1:npar) {
        # put parameters in separate variables
        joe <- paste(pnames[[i]], "<-", pnum[[i]])
        eval(parse(text = joe))
    }
### NEWNLS -- 140716
## ?? Build resfn, jacfn, call nlfb, reporting?
    tresfn<-model2resfun(modelformula, pvec) 
    tjacfn<-model2jacfun(modelformula, pvec) 



    pnum <- as.vector(pnum)
    names(pnum) <- pnames
    result <- list(resid = resbest, jacobian = Jac, feval = feval, 
        jeval = jeval, coefficients = pnum, ssquares = ssbest, lower=lower, upper=upper, maskidx=maskidx)
    class(result) <- "nlmrt"
    result
}
