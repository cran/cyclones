
within <- function(x,y,X,Y,plot=FALSE,test=FALSE) {

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


within.old <- function(x,y,X,Y,test=FALSE) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]   
  n1 <- (x >= min(X,na.rm=TRUE)) & (x <= max(X,na.rm=TRUE)) &
        (y >= min(Y,na.rm=TRUE)) & (y <= max(Y,na.rm=TRUE))
  #print(sum(n1))
  N <- rep(FALSE,length(n1))
  nx <- length(X)

  if (test) points(x[n1],y[n1],pch=20,col="green",cex=0.9)
  first <- TRUE
  #print(sum(n1,na.rm=TRUE))
  if (sum(n1)>0) {
  for (i in seq(1,sum(n1,na.rm=TRUE),by=1)) {
    xx <- x[n1][i]; yy <- y[n1][i]
    if (is.finite(yy) & is.finite(xx)) {
      icross <- c((X[2:nx] - xx)*(X[1:(nx-1)] - xx) <= 0)
      if (sum(icross,na.rm=TRUE)==2) {
        X.x <- X[2:nx][icross]
        Y.x <- Y[2:nx][icross]
        if (length(Y.x)==2) {
          if (yy < min(Y.x,na.rm=TRUE) | yy > max(Y.x,na.rm=TRUE) |
              xx < min(X.x,na.rm=TRUE) | xx > max(X.x,na.rm=TRUE)) {
            N[n1][i] <- FALSE
            if (first) {
              if (test) points(xx,yy,pch=20,col="steelblue",cex=2.5)
              if (test) lines(X.x,Y.x,col="steelblue")
              if (test) lines(c(X.x,xx),c(Y.x,yy),col="steelblue",lty=3)
#              first <- FALSE
            }
          } else {
              if (test) print(paste(i,"point:",xx,yy,"X.x:",X.x,"Y.y",Y.x))
              if (test) points(xx,yy,pch=20,col="yellow",cex=2.5)
              N[n1][i] <- TRUE
          }
        }  else {
        if (test) points(xx,yy,pch=20,col="grey80",cex=2.5)
        n1[n1][i] <- FALSE
        if (test) print(paste(" ---> ",i,"point:",xx,yy,"X.x:",X.x,"Y.y",Y.x))
        } # else & if length==2
      } else {
        if (test) points(xx,yy,pch=20,col="wheat",cex=2.5)
        n1[n1][i] <- FALSE
        if (test) print(paste(" =======> ",i,"point:",xx,yy,"X.x:",X.x,"Y.y",Y.x))
      } # else & if sum==2?
    }  # is.finite? 
  } # for loop i
  } # end if
  invisible(N)
}

test.within <- function(a=3,b=5) {
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


