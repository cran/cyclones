# Calculates the y-derivatives for gridded data in longitude-
# latitude coordinates. 
#
# d/dy = 1/r * d/dPHI 
#
# where PHI is the latitude in radians and THETA the longitude.
# R.E. Benestad, 26.06.2003.
#
# also see reg.cal.R
dY <- function(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE,plot=FALSE,chk.conf=1) {
  ny <- length(lat)
  nx <- length(lon)
  if (is.null(maxhar)) maxhar <- ny
  maxhar <- min(ny,maxhar)
  mask <- !is.finite(Z)
  theta <- pi*lon/180
  phi <- pi*lat/180
  W <- 2*pi/ny
  a <- matrix(rep(0,ny*nx),nx,ny)
  b <- matrix(rep(0,ny*nx),nx,ny)
  c <- matrix(rep(0,ny*nx),nx,ny)
  dZ <- matrix(rep(0,ny*nx),nx,ny)
  Z.fit <- matrix(rep(0,ny*nx),nx,ny)
  dy <- distAB(lon[1],lat[1],lon[1],lat[2])[1]/1000  # REB: fix - units in km;
  attr(dy,"units") <- "km"                        # REB: fix - units in km;
  
  for (i in 1:nx) {
    y <- Z[i,]
    for (iw in 1:maxhar) {
     good <- is.finite(y)
     if (sum(good)>7) {
       wt <- iw*W*seq(1,ny,by=1)
       x1 <- cos(wt); x2 <- sin(wt)
       harmfit <- data.frame(y=y, x1=x1, x2=x2)
       harmonic <- lm(y ~ x1 + x2,data=harmfit)
       y[good] <- harmonic$residual
       c[i,iw] <- harmonic$coefficients[1]; if (!is.finite(c[i,iw])) c[i,iw] <- 0
       a[i,iw] <- harmonic$coefficients[2]; if (!is.finite(a[i,iw])) a[i,iw] <- 0
       b[i,iw] <- harmonic$coefficients[3]; if (!is.finite(b[i,iw])) b[i,iw] <- 0
       if (!is.null(chk.conf)) {
         stats <- summary(harmonic)
         if (abs(stats$coefficients[4])*chk.conf > abs(stats$coefficients[1]))  c[i,iw] <- 0
         if (abs(stats$coefficients[5])*chk.conf > abs(stats$coefficients[2]))  a[i,iw] <- 0
         if (abs(stats$coefficients[6])*chk.conf > abs(stats$coefficients[3]))  b[i,iw] <- 0
       }
       dZ[i,] <- dZ[i,] + iw*W*( -a[i,iw]*sin(wt) + b[i,iw]*cos(wt) ) 
        
#      Z.fit[i,] <- Z.fit[i,] + predict(harmonic,newdata=harmfit)
       Z.fit[i,] <- Z.fit[i,] + a[i,iw]*cos(wt) + b[i,iw]*sin(wt) + c[i,iw]
     }
    }
    if (plot) {plot(Z[i,],main=(paste("dY:",i,"of",nx)),xlab="x",ylab="z"); lines(Z.fit[i,],lwd=2,col="grey")}
  }
  dZ <- dZ/dy
  if (mask.bad) {
    dZ[mask] <- NA
    Z.fit[mask] <- NA
  }
  results <- list(Z=Z,a=a,b=b,c=c,dZ=dZ,Z.fit=Z.fit,lon=lon,lat=lat,dy=dy)
  class(results) <- "map"
  attr(results,"spatial units") <- "km"
  attr(results,"long_name") <- "y-derivative"
  attr(results,"descr") <- "dY()"
  invisible(results)
}
