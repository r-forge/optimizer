#' coef.nlsr: Output model coefficients for nlsr object.
#' 
#' A routine to extract and display the coefficients for a model 
#' estimated by \code{nlxb} or \code{nlfb} in the \code{nlsr} structured
#' \code{object}. 
#' 
#' @section Usage:
#' 
#' coef(object, ...)
#' 
#' \arguments{
#' \item{object}{
#'   An object of class 'nlsr'
#' }
#' \item{\dots}{
#'   Any data needed for the function. We do not know of any!
#' }

#' @param object An object of class 'nlsr'
#' @param ... dot-args to provide exogenous data to the problem
#'     Any data needed for the function. We do not know of any!
#'     
#' @section Details
#' 
#' \code{coef.nlsr} extracts and displays the coefficients for a model 
#'     estimated by \code{nlxb} or \code{nlfb}. 
#'     
#' @section value
#' 
#'  returns the coefficients from the nlsr object. ?? as a named vector
#'  
#' @author John C Nash <nashjc@uottawa.ca>
#' 
#' @seealso 
#'    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
#'    
#' @keyword  nonlinear least squares 
#' 
# coef() function
coef.nlsr <- function(object, ...) {
       out <- object$coefficients
       attr(out,"pkgname")<-"nlsr"
##       invisible(out)
       out # JN 170109
}