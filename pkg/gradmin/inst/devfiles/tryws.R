line <- function(xx, ws){
  ws$tester <- 999
  res <- 2*xx
  res  
}


tryws<-function(x0, control=list(),...) {
cat("control:")
print(control)
npar <- length(x0)
# set up workspace
ctrl <- list(
  maxfevals = npar*500,
  tester = 100,
  trace=0
)  

ncontrol <- names(control)
nctrl <- names(ctrl)
for (onename in ncontrol) {
  if (onename %in% nctrl) {
    ctrl[onename]<-control[onename]
  }
}
cat("tryws::ctrl:")
print(ctrl)
ws <- list2env(ctrl) # Workspace
cat("tryws::ws:")
print(str(ws))
cat("ws$tester =", ws$tester,"\n")
xx <- line(x0,ws)
cat("ws$tester =", ws$tester,"\n")
cat("tryws::xx:")
print(xx)
out <- xx
out 
}

# main
  x00 <- c(4, 3, 2)
  newx <- tryws(x00)
  cat("newx:")
  print(newx)