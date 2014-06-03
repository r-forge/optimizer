# check environments in R
    topenv<-new.env()
top <- function(x, y){
    topenv$g1<-123
    topenv$g2<-321
    topenv$g4<-NA
    joe<-lev2(y, x)
    print(joe)
    print(topenv$g2)
    print(topenv$g1)
}
lev2<-function(yy, xx){
  environment(lev2)<-topenv
  # no topenv here
     jim<-lev3a(xx, yy)
}
lev3a<-function(aa, bb){
    environment(lev2)<-topenv
    big<-aa*topenv$g1+bb*topenv$g2
    topenv$g3<-big*2
    big
}

print(top(1,3))
