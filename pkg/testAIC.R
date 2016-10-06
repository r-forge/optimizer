# testAIC.R
# function to compute AIC etc. for nlmrt

modeldiags <- function(nlsobj) {

  resvec <- nlsobj$resid
  n <- length(resvec)
  par <- nlsobj$coefficients
  k <- length(par) + 1 # Should we add 1, if so why for nonlinear models
  # ?? for linear models this includes the intercept
  rss <- as.numeric(crossprod(resvec))
  # want 2 k - 2 log(L)
  # log(L) = (- n log(2 pi) -n + n log(n) - n log(rss) ) /2
  # 2 * log(L) = n log(2 pi) - n + n log(n) - n log(rss)
  AIC <- 2 * k + n * log(rss) - n * log(n) + n - n * log(2 * pi)  
  BIC <- n * log(rss/n) + k * log(n)
  # others e.g., Cp
  diags <- list(AIC = AIC, BIC = BIC)
}

logLik.nlmrt <- function (object, ...) {
  res <- object$resid
  n <- length(res)
  val <- -n * (log(2 * pi) + 1 - log(n) + log(sum(res^2)))/2
  attr(val, "df") <- 1L + length(coef(object))
  ## JN??: Is this k or degrees of freedom? surely nobs - length(coef(object)), maybe adjusted by 1 ??
  attr(val, "nobs") <- attr(val, "nall") <- n
  class(val) <- "logLik"
  val
}

ex_dat <- structure(list(x = c(60, 45, 25, 20, 12, 7, 4, 1),
   y = c(0.026790358,  0.028404755, 0.028483591, 0.057717406,
   0.073156956, 0.096398739,  0.172890956, 0.635589174)),
   .Names = c("x", "y"), class = "data.frame", row.names = 667:674)

mod1 <- nlxb(y ~ A1*exp(-exp(beta1)*x),
             data = ex_dat,
             start = c(A1 = 1, beta1 = -1))

ll1 <- logLik.nlmrt(mod1)
-2 * as.numeric(ll1) + log(attr(ll1, "nobs")) * attr(ll1, "df")


print(modeldiags(mod1))


mod2 <- nlxb(y ~ A1*exp(-exp(beta1)*x) + A2*exp(-exp(beta2)*x),
             data = ex_dat,
             start = c(A1 = 1, beta1 = -1,
                       A2 = .6, beta2 = -2))

ll2 <- logLik.nlmrt(mod2)
-2 * as.numeric(ll2) + log(attr(ll2, "nobs")) * attr(ll2, "df")

print(modeldiags(mod2))
