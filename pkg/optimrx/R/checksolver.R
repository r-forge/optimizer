checksolver <- function(method){
    basestats <- c("Nelder-Mead","BFGS","L-BFGS-B","CG","SANN", "nlm", "nlminb", "hjn")
    if (method %in% basestats) return(method)

    if (requireNamespace(method, quietly = TRUE)) {
      return(method)
    } else {
      warning("Method ",method," is not available")
      return(NULL)
    }     
    NULL # just in case
}
