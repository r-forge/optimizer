# tDanielWood-nls.R
require(NISTO)
pfname <- "DanielWood"
# try it
counters <- new.env()
test1 <- runoptprob(pfilename=pfname, minmeth="nls")
test1