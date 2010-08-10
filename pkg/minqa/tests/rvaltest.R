library(minqa)

## maxfn <- function(x) 10 - crossprod(x - seq_along(x))^2
minfn <- function(x) crossprod(x - seq_along(x))^2 - 10

x0 <- rep.int(pi, 4)
reschk <- function(res) {
    stopifnot(is.list(res),
              inherits(res, "minqa"),
              names(res) == c("par", "fval", "feval", "ierr", "msg"),
              is.numeric(res$par),
              all.equal(res$par, 1:4, tol = 2e-4),
              is.numeric(res$fval),
              all.equal(as.vector(res$fval), -10, check.attributes = FALSE, tol = 1e-4),
              is.integer(res$feval),
              res$feval > 0)
}
# NOTE: we do not check ierr or msg here. JN 20100810
reschk(ans.nd <- newuoa(x0, minfn, control = list(iprint = 2)))
ans.nd
reschk(ans.ud <- uobyqa(x0, minfn, control = list(iprint = 2)))
ans.ud
reschk(ans.bd <- bobyqa(x0, minfn, control = list(iprint = 2)))
ans.bd

