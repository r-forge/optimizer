library(R6)
library(TeachingDemos)

PolyTrack <- R6Class("PolyTrack",
  public = list(
    parms = list(),
    areas = list(),
    nv = 6,
    PlotIt = TRUE,
    Delay = 0.5,
    nPolys = 5,
    add = function(p,a) {
      i <- length(self$parms) + 1
      self$parms[[i]] <- p
      self$areas[[i]] <- a
      if(self$PlotIt) {
        self$PlotPolys()
        Sys.sleep(self$Delay)
      }
      return(a)
    },
    PlotPolys = function(i=-1) {
      if(i<0) i <- length(self$parms)
      if(i==0) return()
      cols <- hsv(0.6, (1:self$nPolys)/self$nPolys, 1)
      start <- pmax(1, i-self$nPolys+1)
      plotParms <- self$parms[seq(start,i)]
      n <- length(plotParms)
      if(n < self$nPolys) cols <- tail(cols, n)
      coords <- lapply(plotParms, function(x) polypar2XY(self$nv, x))
      plot.new()
      plot.window( xlim=do.call(range, lapply(coords, function(xy) xy$x)),
                   ylim=do.call(range, lapply(coords, function(xy) xy$y)),
                   asp=1)
      for(i in seq_len(n)) {
        polygon(coords[[i]]$x, coords[[i]]$y, border=cols[i], lwd=3)
      }
    }
  )
                     
)


pt1 <- PolyTrack$new()


polyobj <- function(x, penfactor=0) {
  # negative area + penfactor*(sum(squared violations))
  nv = (length(x)+3)/2 # number of vertices
  f <-  -polyarea(nv, x) # negative area
  XY <- polypar2XY(nv, x)
  dist2 <- polydistXY(nv, XY)
  viol <- dist2[which(dist2 > 1)] - 1.0
  f <- f + penfactor * sum(viol)
  pt1$add(x,f)
}

library(minqa)
ub <- c(rep(1,(nv-1)), rep(0.75*pi, (nv-2))) # approx for angles
lb <- c(rep(0, (2*nv-3)))

start<-myhex$par

sol <- bobyqa(start, polyobj, lower=lb, upper=ub, control=list(iprint=3), penfactor=10)



# Redo the plots/animation after the optimization
tkexamp(pt1$PlotPolys(), list(i=list('animate', init=1, from=1, to=length(pt1$parms), delay=pt1$Delay*100)))



pt2 <- PolyTrack$new()

polyobj <- function(x, penfactor=0) {
  # negative area + penfactor*(sum(squared violations))
  nv = (length(x)+3)/2 # number of vertices
  f <-  -polyarea(nv, x) # negative area
  XY <- polypar2XY(nv, x)
  dist2 <- polydistXY(nv, XY)
  viol <- dist2[which(dist2 > 1)] - 1.0
  f <- f + penfactor * sum(viol)
  pt2$add(x,f)
}


library(Rvmmin)
solvm<-Rvmmin(start, polyobj, lower=lb, upper=ub, control=list(trace=1), penfactor=10)

# Redo the plots/animation after the optimization
tkexamp(pt2$PlotPolys(), list(i=list('animate', init=1, from=1, to=length(pt2$parms), delay=pt2$Delay*100)))
