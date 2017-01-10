## weighted nonlinear regression
Treated <- Puromycin[Puromycin$state == "treated", ]
weighted.MM <- function(resp, conc, Vm, K)
{
  ## Purpose: exactly as white book p. 451 -- RHS for nls()
  ##  Weighted version of Michaelis-Menten model
  ## ----------------------------------------------------------
  ## Arguments: 'y', 'x' and the two parameters (see book)
  ## ----------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Mar 2001
  
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}

Pur.wtpred <- nls( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
               start = list(Vm = 200, K = 0.1))
summary(Pur.wtpred)

Pur.wtw <- nls( rate ~ (Vm * conc)/(K + conc), data = Treated,
                  start = list(Vm = 200, K = 0.1), weights=1/rate)
summary(Pur.wtw)
library(nlsr)
Pur.wtx <- nlxb( rate ~ (Vm * conc)/(K + conc), data = Treated,
                start = list(Vm = 200, K = 0.1), weights=1/Treated$rate)
summary(Pur.wtx)
stw <- c(Vm = 200, K = 0.1)
wres.MM <- function(prm, rate, conc) {
      print(prm)
      print(rate)
      print(conc)
      Vm <- prm[1]
      K <- prm[2]
      pred <- (Vm * conc)/(K + conc)
    (rate - pred) / sqrt(rate) # NOTE: NOT pred here
}

# test
print(wres.MM(stw, rate=Treated$rate, conc=Treated$conc))

Pur.wtf <- nlfb(start=stw,  resfn=wres.MM, rate = Treated$rate, conc=Treated$conc, 
                  trace=TRUE) # Note: no weights needed here
summary(Pur.wtf)
print(Pur.wtf)

wss.MM <- function(prm, rate, conc) {
  ss <- as.numeric(crossprod(wres.MM(prm, rate, conc)))
}
library(optimrx)
osol <- optimr(stw, fn=wss.MM, gr="grnd", method="Nelder-Mead", control=list(trace=1), 
            rate=Treated$rate, conc=Treated$conc)
print(osol)
