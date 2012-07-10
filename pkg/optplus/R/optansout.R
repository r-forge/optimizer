optansout <- function(ansdf, filename) {
    ##### OPEN ISSUES: (any date order)
    # 120709 Not printing parameters == but does when source'd!    
    ##### IMPLEMENTED: (reverse date order)
    
    # A funtion to display and print to file (if present) the
    #   output of optimx
    if (!is.null(filename)) {
        sink(filename)
    }
    # if (! exists(filename)) { sink(filename) } # may need to
    #   check paths
    tpar <- ansdf$par
    tdf <- ansdf
#    ltvec <- length(tpar[[1]]) 
    for (i in 1:length(tpar)) {
        tvec <- tpar[[i]]
        ltvec <- length(tvec)
        if (ltvec > 5) {
            tvec <- tvec[1:5]
            cat("truncating parameters\n")
        } else cat("no shortening of parameters\n")
        tpar[[i]] <- tvec
    }
    tdf$par <- tpar
    if (ltvec > 5) {
        names(tdf)[2] <- "first.5.par"
    } else {
        names(tdf)[2] <- "par"
    }
    print(tdf)
    if (!is.null(filename)) {
        sink()
    }
    rm(tdf)
    rm(tpar)
    return(TRUE)
}
