##------------------------------------------------------------------
require(circular) ## for Bessel function I.0

## Data:
dd <- c(0.9975948929787, 0.9093316197395, 0.7838819026947, 
0.9096108675003, 0.8901804089546, 0.2995955049992, 0.9461286067963, 
0.8248071670532, 0.2442084848881, 0.2836948633194, 0.7353935241699, 
0.5812761187553, 0.8705610632896, 0.8744471669197, 0.7490273118019, 
0.9947383403778, 0.9154829382896, 0.8659985661507, 0.6448246836662, 
0.8588128685951, 0.7347437739372, -0.1645197421312, 0.970999121666, 
0.8038327097893, 0.9558997154236, 0.6846113204956, 0.6286814808846, 
0.9201356172562, 0.9422197341919, 0.3470877110958, 0.4154576957226, 
0.0721184238791, 0.14151956141, -0.6142936348915, -0.4688512086868, 
0.6805665493011, 0.3594025671482, 0.8991097211838, 0.7656877636909, 
0.9282909035683, 0.9454715847969, 0.9766132831573, 0.4316343963146, 
0.62679708004,   0.2093886137009, 0.3937581181526, 0.4254160523415, 
0.8684504628181, 0.3844584524632, 0.9578431844711, 0.956972181797, 
0.4456568360329, 0.9793710708618, 0.5825698971748, 0.929228246212, 
0.9211971759796, 0.9407976865768, 0.821156680584, 0.2048042863607, 
0.6473184227943, 0.9456319212914, 0.7021154165268, 0.9761978387833, 
0.1485801786184, 0.2195029109716, 0.5378785729408, 0.8304615020752, 
0.8596342802048, 0.950027525425, 0.9102076888084, 0.5108731985092, 
0.7200184464455, 0.3571084141731, 0.9765330553055, -0.143017962575, 
0.8576183915138, 0.1283493340015, -0.3226098418236, 0.7031792402267, 
0.8708582520485, 0.56754809618, 0.060470353812, 0.8015220761299, 
0.7363410592079, 0.671902179718, 0.8082517385483, 0.9468197822571, 
0.9729647636414, 0.7919752597809, 0.9539568424225, 0.4840737581253, 
0.850653231144, 0.5909016132355, 0.8414449691772, 0.9699150323868)

xlims <- c(-1,1)
bw <- 0.05
b <- seq(xlims[1],xlims[2],by=bw)   ;  nb <- length(b)
h <- hist( dd, breaks=b, plot=FALSE)

FisherAvgdPdf <- function(theta,theta0,kappa){
   A <- kappa/(2*sinh(kappa))
   A * I.0( kappa*sin(theta)*sin(theta0) ) * exp( 
kappa*cos(theta)*cos(theta0) )
}

# JN added ans<- and trace=TRUE rather than warnonly=TRUE
ans<-nls(dens ~ FisherAvgdPdf(theta,theta0,kappa),
     data = data.frame( theta=acos(h$mids), dens=h$density ),
     start=c( theta0=0.5, kappa=4.0 ),
     algorithm="port", lower=c(0.0001,0.0001), upper=c(acos(xlims[1]),500),
     control=list(trace=TRUE) )

print(ans)
require("optimx")
xx<-c(0.5,4)
theta<-acos(h$mids)
dens<-h$density
print(theta)
print(dens)
df<-data.frame(dens, theta)

#resbett<-function(xx){
resbett<-function(xx, df){
# residual function for Bett problem at data theta)
   t0<-xx[1]
   kk<-xx[2]
   dens<-df$dens
   theta<-df$theta
   r <- dens - FisherAvgdPdf(theta, t0, kk)
}

# fbett<-function(xx, theta, dens){
fbett<-function(xx, df){
    res<-resbett(xx, df)
    f<-as.numeric(crossprod(res))
}
res0<-resbett(xx, df)
cat("Initial res:")
print(res0)
cat("initial fn:")
f0<-fbett(xx, df)
print(f0)


mymeth<-c("Nelder-Mead", "BFGS", "L-BFGS-B", "Rvmmin", "Rcgmin", "spg", 
             "ucminf", "nlminb", "nlm")

aopt<-optimx(xx, fn=fbett, lower=c(0.0001,0.0001), upper=c(acos(xlims[1]),500),
               method=mymeth, control=list(trace=1), df=df)
#aopt<-optimx(xx, fn=fbett, control=list(trace=1), df=df)
print(aopt)
   

sessionInfo()
##------------------------------------------------------------------

