model2jacfun <- function(resformula, pvec, filename=NULL) {
   pnames<-names(pvec)
# Creates Jacobian function
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
#   cat("xx2:")
#   print(xx2)
#   cat("xxparm:")
#   print(xxparm)
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
   cat("vstr:")
   print(vstr)
   tmp<-readline("...")
    ff <- vector("list", length(xx2))
   names(ff) <- xx2
   parts<-strsplit(as.character(es), "~")[[1]]
   if (length(parts)!=2) stop("Model expression is incorrect!")
   lhs<-parts[1]
   rhs<-parts[2]
#  And build the residual at the parameters
   resexp<-paste(rhs,"-",lhs, collapse=" ") # build the residuals
   jacexp<-deriv(parse(text=resexp), pnames) # gradient expression
#   cat("jacexp:")
#   print(jacexp)
#   tmp<-readline("cont.")
   dvstr<-""
   if (nvar>0) {
      for (i in 1:nvar){
         dvstr<-paste(dvstr,xxvars[i], sep='')
         if (i<nvar) dvstr<-paste(dvstr,", ",sep='')
      }
   }
   jfstr<-paste("localdf<-data.frame(",dvstr,");\n", sep='')
   jfstr<-paste(jfstr,"jstruc<-with(localdf,eval(",jacexp,"))", sep="") ##3
   pparse<-""
   for (i in 1:npar){
      pparse<-paste(pparse, "   ",pnames[[i]],"<-prm[[",i,"]]\n", sep='')
   }
#   myjstr<-paste("myjac<-function(prm, ",vstr,"){\n",
   myjstr<-paste("{\n",
      pparse,"\n ",
      jfstr," \n",
      "jacmat<-attr(jstruc,'gradient')\n ",
      "return(jacmat)\n }",sep='')
   if (! is.null(filename)) write(myjstr, file=filename) # write out the file
   tparse<-try(body(myjac)<-parse(text=myjstr)) 
   # This may be cause trouble if there are errors
   if (class(tparse) == "try-error") stop("Error in Jacobian code string")
   return(myjac)      
}


