# Calculates the x-derivatives for gridded data in longitude-
# latitude coordinates. After Gill (1982) p. 94
#
# dA/dx = 1/(r cos(PHI)) * d/d THETA  
#
# where PHI is the latitude in radians and THETA the longitude.
# R.E. Benestad, 26.06.2003.
#
# also see reg.cal.R

dX <- function(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE) {
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
  dZ <- matrix(rep(0,ny*nx),nx,ny)
  Z.fit <- matrix(rep(0,ny*nx),nx,ny)
  for (j in 1:ny) {
    y <- Z[,j]
    for (iw in 1:maxhar) {
     good <- is.finite(y)
     if (sum(good)>10) {
      wt <- iw*W*seq(1,nx,by=1)
      x1 <- cos(wt); x2 <- sin(wt)
      harmfit <- data.frame(y=y, x1=x1, x2=x2)
      harmonic <- lm(y ~ x1 + x2,data=harmfit)
      y[good] <- harmonic$residual
      c[j,iw] <- harmonic$coefficients[1]
      a[j,iw] <- harmonic$coefficients[2]
      b[j,iw] <- harmonic$coefficients[3]
#      Z.fit[,j] <- Z.fit[,j] + predict(harmonic,newdata=harmfit)
      Z.fit[,j] <- Z.fit[,j] + a[j,iw]*cos(wt) + b[j,iw]*sin(wt) + c[j,iw]
      dZ[,j] <- dZ[,j] +iw*W*( -a[j,iw]*sin(wt) + b[j,iw]*cos(wt) )
    }
   }
#    print(c(r,cos(phi[j]),j,r*cos(phi[j])))
    dZ[,j] <- dZ[,j]/(r*cos(phi[j]))
  }
  if (mask.bad) {
    dZ[mask] <- NA
    Z.fit[mask] <- NA
  }
  results <- list(Z=Z,a=a,b=b,c=c,dZ=dZ,Z.fit=Z.fit,lon=lon,lat=lat)
  class(results) <- "map"
  attr(results,"long_name") <- "y-derivative"
  attr(results,"descr") <- "dX.R"
  invisible(results)
}
