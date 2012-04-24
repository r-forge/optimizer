model2resfun <- function(resformula, pvec, filename) {
   xx <- all.vars(resformula)
   pnames<-names(pvec)
   if (is.null(pnames) ) stop("MUST have named parameters in pvec")
   rp <- match(pnames, xx) # Match names to parameters
   # ?? How to ensure there are names?
   xx2 <- c(xx[rp], xx[-rp])
   xxparm<-xx[rp]
   pstr<-"c("
   npar<-length(xxparm)
   if (npar>0) {
      for (i in 1:npar){
         pstr<-paste(pstr,"\"",xxparm[i],"\"", sep='')
         if (i<npar) pstr<-paste(pstr,", ",sep='')
      }
   }
   pstr<-paste(pstr,")",sep='')
   xxvars<-xx[-rp]
   nvar<-length(xxvars)
   vstr<-""
   if (nvar>0) {
      for (i in 1:nvar){
         vstr<-paste(vstr,xxvars[i]," = NULL", sep='')
         if (i<nvar) vstr<-paste(vstr,", ",sep='')
      }
   }
   ff <- vector("list", length(xx2))
   names(ff) <- xx2
   sf<-as.character(resformula)
   if ((length(sf)!=3) && (sf[1]!="~")) stop("Bad model formula expression")
   lhs<-sf[2] # NOTE ORDER formula with ~ puts ~, lhs, rhs
   rhs<-sf[3]
   resexp<-paste(rhs,"-",lhs, collapse=" ") # build the residuals
   fnexp<-paste("eval(crossprod(",resexp,"))", sep="") ##3
   pparse<-paste("for (i in 1:length(prm) ){\n",
      "joe<-paste(names(prm)[[i]],\"<-\",prm[[i]]);\n",
      " eval(parse(text=joe));\n",
      "}", sep='')
   myfstr<-paste("myfn<-function(prm, ",vstr,"){;\n",
      "if ( is.null(names(prm)) ) { warning('Parameter vector must have names'); return(NA)};\n",
      pparse,";\n ",
      fnexp,"\n }",sep='')
   write(myfstr, file=filename) # write out the file
   myfn<-source(file=filename)$value # This may be inefficient, but ...
   attr(myfn, "source-text")<-myfstr
   return(myfn)      
}


