## Put in R directory. 
if(getRversion() >= "2.15.1") { utils::globalVariables(c('envroot')) } # Try declaring here 
groot<-list(ifn=0, igr=0, ftrace=FALSE, fn=NA, gr=NA, label="none")
envroot <- list2env(groot) # Note globals in FnTrace
