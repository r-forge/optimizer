
 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()

#require("stats")  

fuzz <- 1e-6 
all.ok <- TRUE  


# test MASS 4th ed. p 322-324

  data(ability.cov)
  ability.FA <- factanal(factors = 1, covmat=ability.cov)

# check if this was suppose to do something in MASS
#  loadings(rotation(ability.FA, rotation="oblimin"))


  ability.FA2 <- factanal(factors = 2, covmat = ability.cov)
# check if this was suppose to do something in MASS
#  loadings(rotation(ability.FA2, rotation="oblimin"))

  oblimin(loadings(ability.FA2))

# help for factanal has
# example(ability.cov)

# factanal(factors = 2, covmat = ability.cov, rotation="oblimin")
# factanal(factors = 2, covmat = ability.cov, rotation="oblimin",
#          control=list(opt=))

# promax(ability.FA2$loadings)

# varimax(ability.FA2$loadings)

# factanal(factors = 2, covmat = ability.cov, scores = Bartlett, rotation="oblimin")


# tst <- t(matrix(c(), 6, 2))
# 
# if( fuzz < max(abs(... - tst ))) {
#    cat("Calculated value is not the same as test value in test x. Value:\n")
#    print(..., digits=18)
#    cat("difference:\n")
#    print(... - tst, digits=18)
#    all.ok <- FALSE  
#    } 

# from help for Harman74.cor
#     Harman, H. H. (1976) _Modern Factor Analysis_, Third Edition
#     Revised, University of Chicago Press, Table 7.4.
 
     require(stats)
 
     (Harman74.FA <- factanal(factors = 1, covmat = Harman74.cor))
     for(factors in 2:5) print(update(Harman74.FA, factors = factors))
     Harman74.FA <- factanal(factors = 5, covmat = Harman74.cor,
                             rotation="promax")
     print(Harman74.FA$loadings, sort = TRUE)
 

cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
