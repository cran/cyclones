# Calculates the y-derivatives for gridded data in longitude-
# latitude coordinates. 
#
# d/dy = 1/r * d/dPHI 
#
# where PHI is the latitude in radians and THETA the longitude.
# R.E. Benestad, 26.06.2003.
#
# also see reg.cal.R
dY <- function(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE) {
  ny <- length(lat)
  nx <- length(lon)
  if (is.null(maxhar)) maxhar <- ny
  maxhar <- min(ny,maxhar)
  mask <- !is.finite(Z)
  theta <- pi*lon/180
  phi <- pi*lat/180
  W <- pi/ny
  a <- matrix(rep(0,ny*nx),nx,ny)
  b <- matrix(rep(0,ny*nx),nx,ny)
  c <- matrix(rep(0,ny*nx),nx,ny)
  dZ <- matrix(rep(0,ny*nx),nx,ny)
  Z.fit <- matrix(rep(0,ny*nx),nx,ny)
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
       c[i,iw] <- harmonic$coefficients[1]
       a[i,iw] <- harmonic$coefficients[2]
       b[i,iw] <- harmonic$coefficients[3]
       dZ[i,] <- dZ[i,] + iw*W*( -a[i,iw]*sin(wt) + b[i,iw]*cos(wt) ) 
        
#      Z.fit[i,] <- Z.fit[i,] + predict(harmonic,newdata=harmfit)
       Z.fit[i,] <- Z.fit[i,] + a[i,iw]*cos(wt) + b[i,iw]*sin(wt) + c[i,iw]
     }
    }
    dZ[i,] <- dZ[i,]/r
  }
  if (mask.bad) {
    dZ[mask] <- NA
    Z.fit[mask] <- NA
  }
  results <- list(Z=Z,a=a,b=b,c=c,dZ=dZ,Z.fit=Z.fit,lon=lon,lat=lat)
  class(results) <- "map"
  attr(results,"long_name") <- "y-derivative"
  attr(results,"descr") <- "dY.R"
  invisible(results)
}
