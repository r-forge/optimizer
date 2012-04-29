model2resfun <- function(resformula, pvec, filename=NULL) {
   pnames<-names(pvec)
#   cat("pnames:")
#   print(pnames)
   if (is.null(pnames) ) stop("MUST have named parameters in pvec")
   if (is.character(resformula)){
      es<-resformula
   } else {
      tstr<-as.character(resformula) # note ordering of terms!
      es<-paste(tstr[[2]],"~",tstr[[3]],'')
   }
   xx <- all.vars(parse(text=es))
#   cat("xx:")
#   print(xx)
   rp <- match(pnames, xx) # Match names to parameters
# ?? How to ensure there are names?
   xx2 <- c(xx[rp], xx[-rp])
   xxparm<-xx[rp]
   cat("xx2:")
   print(xx2)
   cat("xxparm:")
   print(xxparm)
   pstr<-"c("
   npar<-length(xxparm)
   if (npar>0) {
      for (i in 1:npar){
         pstr<-paste(pstr,'"',xxparm[i],'"', sep='')
         if (i<npar) pstr<-paste(pstr,", ",sep='')
      }
   }
   pstr<-paste(pstr,")",sep='')
   cat("pstr:")
   print(pstr)
   tmp<-readline("...")
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
   parts<-strsplit(as.character(es), "~")[[1]]
   if (length(parts)!=2) stop("Model expression is incorrect!")
   lhs<-parts[1]
   rhs<-parts[2]
#  And build the residual at the parameters
   resexp<-paste(rhs,"-",lhs, collapse=" ") # build the residuals
   fnexp<-paste("resids<-as.numeric(eval(",resexp,"))", sep="") ##3
   pparse<-""
   for (i in 1:npar){
      pparse<-paste(pparse, "   ",pnames[[i]],"<-prm[[",i,"]]\n", sep='')
# for diagnostic
#      pparse<-paste(pparse, "   cat(\"",pnames[[i]],"=\",",pnames[[i]],",\"\n\")\n", sep='')
   }
# for diagnostic
#   pparse<-paste(pparse, "cat(\"tt:\")\n  print(tt)\n  cat(\"y:\")\n print(y)\n", sep='')
#   myfstr<-paste("myres<-function(prm, ",vstr,"){\n",
    myfstr<-paste("{\n",
      pparse,"\n ",
      fnexp,"\n }",sep='')
   if (! is.null(filename)) write(myfstr, file=filename) # write out the file
   tparse<-try(body(myres)<-parse(text=myfstr)) 
   # This may be cause trouble if there are errors
   if (class(tparse) == "try-error") stop("Error in residual code string")
   print(myres)
   return(myres)      
}

