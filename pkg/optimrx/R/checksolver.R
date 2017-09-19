checksolver <- function(method, allmeth, allpkg){
#    basestats <- c("Nelder-Mead","BFGS","L-BFGS-B","CG","SANN", "nlm", "nlminb", "hjn")
#    if (method %in% basestats) return(method)

    imeth <- which(method == allmeth)
    pkg <- allpkg[imeth]

    if (requireNamespace(pkg, quietly = TRUE)) {
      return(method)
    } else {
      warning("Package ",pkg," for method ",method," is not available")
      return(NULL)
    }     
    NULL # just in case
}
