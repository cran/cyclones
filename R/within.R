
within <- function(x,y,X,Y,test=FALSE,plot=FALSE) {

  center.X <- mean(X,na.rm =TRUE)
  center.Y <- mean(Y,na.rm =TRUE)
  
  angle.xy <- atan2(x-center.X,y-center.Y)
  dist.xy <- sqrt( (x-center.X)^2 + (y-center.Y)^2 )

  angle.XY <- atan2(X-center.X,Y-center.Y)  
  dist.XY <- sqrt( (X-center.X)^2 + (Y-center.Y)^2 )

  dist.hat <- approx(angle.XY,dist.XY,angle.xy,rule=2)$y
  within <- ( dist.xy <= dist.hat )

  within[is.na(within)] <- FALSE

  if (plot) {
   plot(angle.XY,dist.XY,type="l")
   points(angle.xy,dist.xy,col="blue")
   lines(c(t(cbind(angle.xy[within],angle.xy[within],rep(NA,sum(within))))),
         c(t(cbind(dist.hat[within],dist.xy[within], rep(NA,sum(within))))))
   points(angle.xy[within],dist.hat[within],pch=20)
   points(angle.xy[within],dist.xy[within],pch=20)
  }
  #points(x[within],y[within],cex=4,col="blue")
  within 
}



testWithin <- function(a=3,b=5) {
  x <- seq(-a,a,length=1000)
  y <- b*sqrt(1 - (x/a)^2)
  plot(c(x,-x),c(-y,y),type="l",xlim=c(-a-b,a+b),ylim=c(-a-b,a+b),
       main="test.within")
 if (b>a) {cosT <- cos(atan2(a,b)); sinT <- sin(atan2(a,b))} else
          {cosT <- cos(atan2(b,a)); sinT <- sin(atan2(b,a))}
  R <- matrix(c(cosT,sinT,-1*sinT,cosT),2,2)
  reg1 <- cbind(x,y); reg1 <- R %*% t(reg1)
  x.reg1 <- reg1[1,];  y.reg1 <- reg1[2,]
  lines(c(-x.reg1,x.reg1),c(y.reg1,-y.reg1),lwd=3,col="grey70")
  z1 <- rnorm(1000)*a; z2 <- rnorm(1000)*b
  ii <- within(z1,z2,c(-x.reg1,x.reg1),c(y.reg1,-y.reg1),test=TRUE)
  points(z1,z2,pch=20,col="grey20",cex=0.7)
  points(z1[ii],z2[ii],pch=21,col="red")

  lines(c(-a,a)*3,c(-b,b)*3,lwd=2,lty=2,col="green")
  lines(c(a,-a)*3,c(-b,b)*3,lwd=2,lty=2,col="darkblue")
}


