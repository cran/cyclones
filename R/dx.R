# Calculates the x-derivatives for gridded data in longitude-
# latitude coordinates. After Gill (1982) p. 94
#
# dA/dx = 1/(r cos(PHI)) * d/d THETA  
#
# where PHI is the latitude in radians and THETA the longitude.
# R.E. Benestad, 26.06.2003.
#
# also see reg.cal.R

dX <- function(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE,plot=FALSE,chk.conf=1,accuracy=NULL) {
# REB: 08.02.2007: To improve accuracy/spatial resolution of the fit
  if (is.null(accuracy)) accuracy <-  max(diff(lon))
  LON <- seq(min(lon),max(lon),by=accuracy)
  NX <- length(LON)

  ny <- length(lat)
  nx <- length(lon)
  
  if (is.null(maxhar)) maxhar <- nx
  maxhar <- min(nx,maxhar)
  theta <- pi*lon/180
  phi <- pi*lat/180
  mask <- !is.finite(Z)
  W <- 2*pi/nx
  a <- matrix(rep(0,ny*nx),ny,nx)
  b <- matrix(rep(0,ny*nx),ny,nx)
  c <- matrix(rep(0,ny*nx),ny,nx)
#  dZ <- matrix(rep(0,ny*nx),nx,ny)
#  Z.fit <- matrix(rep(0,ny*nx),nx,ny)
# REB: 08.02.2007:
  dZ <- matrix(rep(0,ny*NX),NX,ny)
  Z.fit <- matrix(rep(0,ny*NX),NX,ny)
  dx <- rep(NA,ny); attr(dx,"units") <- "km"
  
  for (j in 1:ny) {
    y <- Z[,j]
    dx[j] <- distAB(lon[1],lat[j],lon[2],lat[j])/1000  # REB: fix - units in km
    for (iw in 1:maxhar) {
     good <- is.finite(y)
     if (sum(good)>10) {
      wt <- iw*W*seq(1,nx,by=1)
      WT <- seq(min(wt),max(wt),by=accuracy)    # REB: 08.02.2007
      x1 <- cos(wt); x2 <- sin(wt)
      harmfit <- data.frame(y=y, x1=x1, x2=x2)
      harmonic <- lm(y ~ x1 + x2,data=harmfit)
      y[good] <- harmonic$residual
      c[j,iw] <- harmonic$coefficients[1]; if (!is.finite(c[j,iw])) c[j,iw] <- 0
      a[j,iw] <- harmonic$coefficients[2]; if (!is.finite(a[j,iw])) a[j,iw] <- 0
      b[j,iw] <- harmonic$coefficients[3]; if (!is.finite(b[j,iw])) b[j,iw] <- 0
      if (!is.null(chk.conf)) {
         stats <- summary(harmonic)
         if (abs(stats$coefficients[4])*chk.conf > abs(stats$coefficients[1]))  c[j,iw] <- 0
         if (abs(stats$coefficients[5])*chk.conf > abs(stats$coefficients[2]))  a[j,iw] <- 0
         if (abs(stats$coefficients[6])*chk.conf > abs(stats$coefficients[3]))  b[j,iw] <- 0
      }
      Z.fit[,j] <- Z.fit[,j] + a[j,iw]*cos(wt) + b[j,iw]*sin(wt) + c[j,iw]
      dZ[,j] <- dZ[,j] +iw*W*( -a[j,iw]*sin(wt) + b[j,iw]*cos(wt) )
# REB: 08.02.2007:
#      Z.fit[,j] <- Z.fit[,j] + a[j,iw]*cos(WT) + b[j,iw]*sin(WT) + c[j,iw]
#      dZ[,j] <- dZ[,j] +iw*W*( -a[j,iw]*sin(WT) + b[j,iw]*cos(WT) ) 

    }
   }
    if (plot) {plot(Z[,j],main=paste("dX:",j,"of",ny),xlab="y",ylab="z")
               lines(Z.fit[,j],lwd=2,col="grey")}
#    dZ[,j] <- dZ[,j]/(r*cos(phi[j]))
    dZ[,j] <- dZ[,j]/dx[j]
  }
  if (mask.bad) {
    dZ[mask] <- NA
    Z.fit[mask] <- NA
  }
  results <- list(Z=Z,a=a,b=b,c=c,dZ=dZ,Z.fit=Z.fit,lon=LON,lat=lat,dx=dx,span=range(lon))
  class(results) <- "map"
  attr(results,"long_name") <- "y-derivative"
  attr(results,"spatial units") <- "km"
  attr(results,"descr") <- "dX()"
  invisible(results)
}
