library(clim.pact)

test.dX <- function(field=NULL,maxhar=7) {
  if (is.null(field)) data(DNMI.slp,envir=environment()) else {
                      DNMI.slp <- field; rm(field) }
  
  slpmap <- t(DNMI.slp$dat[length(DNMI.slp$tim),,])
  lon <- DNMI.slp$lon; lat <- DNMI.slp$lat
  dy.dx <- dX(lon,lat,slpmap,maxhar=maxhar)

  ii <- (1:length(lat))[lat >= mean(lat)][1]; nx <- length(lon)
  Y <- c(DNMI.slp$dat[,ii,])
  #print(c(ii,nx,length(lon)))
  
  newFig()
  x<- dy.dx$dZ[,ii]; y <- diff(dy.dx$Z.fit[,ii])/dy.dx$dx[ii]
  print(c(length(x),length(y), NA, nx,length(lon), NA, ii,length(lat),dy.dx$dx[ii]))
  x <- 0.5*(x[1:(nx-1)] + x[2:nx])
  print(summary(x)); print(summary(y))
  plot(x,y,pch=20,col="grey20",xlab="dX()",ylab="diff(best-fit)/dx",
       main="test.dX: dX() versus diff() * dx")

  fit <- lm(y ~ x)
  abline(fit,col="red",lty=2)
  X <- distAB(lon[1],lat[ii],lon,lat[ii])/1000; attr(x,'units') <- 'km'
#  X <- seq(0,1,length=nx)  

  newFig()
  plot(X,slpmap[,ii],pch=20,col="grey20",
       main="test.dX: Reconstruction based on dX() and diff(best-fit)*dx")
  lines(X,dy.dx$Z.fit[,ii],col="red",lwd=2)
  y <- dy.dx$dZ[,ii]*dy.dx$dx[ii]
  y <- dy.dx$Z.fit[1,ii] + (cumsum(y)-y[1])
  scaling <- sd(dy.dx$Z.fit[,ii],na.rm=TRUE)/sd(y,na.rm=TRUE)
  my <- y[1]
  print(paste("sd difference is:",scaling))
  # TEST y <- (y - my)* scaling + my
  lines(X,y,col="darkblue",lty=2,lwd=2)

  Y <- c(Y[1],Y[1]+cumsum(diff(dy.dx$Z.fit[,ii])))
  #print(Y)
  lines(X,Y,col="grey",lty=2)

  YY <- c(Y[1],Y[1]+cumsum(diff(predict(fit)*dy.dx$dx[ii])))
  lines(X[1:length(YY)],predict(fit),col="green",lty=2)
  dy.dx$lm.fit <- fit
  print(summary(fit))
  invisible(dy.dx)
}


test.dY <- function(field=NULL,maxhar=5) {
  if (is.null(field)) data(DNMI.slp,envir=environment()) else {
                      DNMI.slp <- field; rm(field) }
  
  slpmap <- t(DNMI.slp$dat[length(DNMI.slp$tim),,])
  lon <- DNMI.slp$lon; lat <- DNMI.slp$lat
  dy.dx <- dY(lon,lat,slpmap,maxhar=maxhar)
    
  ii <- (1:length(lon))[lon >= mean(lon)][1]; ny <- length(lat)
  Y <- c(DNMI.slp$dat[,,ii])
  #print(c(ii,nx,length(lon)))
  
  newFig()
  x<- dy.dx$dZ[ii,]; y <- diff(dy.dx$Z.fit[ii,])/dy.dx$dy
  print(c(length(x),length(y), NA, ny,length(lat), NA, ii,dy.dx$dy))
  x <- 0.5*(x[1:(ny-1)] + x[2:ny])
  print(summary(x)); print(summary(y))
  plot(x,y,pch=20,col="grey20",xlab="dY()",ylab="diff(best-fit)/dy",
       main="test.dX: dY() versus diff() * dy")

  fit <- lm(y ~ x)
  abline(fit,col="red",lty=2)
  X <- distAB(lon[ii],lat[1],lon[ii],lat)/1000; attr(x,'units') <- 'km'
#  X <- seq(0,1,length=ny)  

  newFig()
  plot(X,slpmap[ii,],pch=20,col="grey20",
       main="test.dY: Reconstruction based on dY() and diff(best-fit)*dy")
  lines(X,dy.dx$Z.fit[ii,],col="red",lwd=2)
  y <- dy.dx$dZ[ii,]*dy.dx$dy
  y <- dy.dx$Z.fit[ii,1] + (cumsum(y)-y[1])
  scaling <- sd(dy.dx$Z.fit[ii,],na.rm=TRUE)/sd(y,na.rm=TRUE)
  my <- y[1]
  print(paste("sd difference is:",scaling))
  # TEST y <- (y - my)* scaling + my
  lines(X,y,col="darkblue",lty=2,lwd=2)

  Y <- c(Y[1],Y[1]+cumsum(diff(dy.dx$Z.fit[ii,])))
  #print(Y)
  lines(X,Y,col="grey",lty=2)

  YY <- c(Y[1],Y[1]+cumsum(diff(predict(fit)*dy.dx$dy)))
  lines(X[1:length(YY)],predict(fit),col="green",lty=2)
  dy.dx$lm.fit <- fit
  print(summary(fit))
  invisible(dy.dx)
}


test.dT <- function(y=NULL,maxhar=15) {
  if (is.null(y)) { data(oslo.t2m,envir=environment()); y <- c(t(oslo.t2m$val)) }
  nt <- length(y)
  dy.dx <- dT(y,maxhar=maxhar)

  newFig()
  x<- dy.dx$dy; y <- diff(dy.dx$y.fit)
  print(c(length(x),length(y), NA, nt))
  x <- 0.5*(x[1:(nt-1)] + x[2:nt])
  print(summary(x)); print(summary(y))
  plot(x,y,pch=20,col="grey20",xlab="dX()",ylab="diff(best-fit)/dx",
       main="test.dT: dT() versus diff() * dt")
  grid()
  fit <- lm(y ~ x)
  abline(fit,col="red",lty=2)

  X <- seq(0,1,length=nt); dx <- 1/nt  
  newFig()
  plot(X,y,pch=20,col="grey20",
       main="test.dX: Reconstruction based on dX() and diff(best-fit)*dx")
  lines(X,dy.dx$y.fit,col="red",lwd=2)
  y.t <- dy.dx$dy*dx
  y.t <- dy.dx$y.fit[1] + (cumsum(y.t)-y.t[1])
  scaling <- sd(dy.dx$y.fit,na.rm=TRUE)/sd(y.t,na.rm=TRUE)
  my <- y.t[1]
  print(paste("sd difference is:",scaling))
  # TEST y <- (y - my)* scaling + my
  lines(X,y.t,col="darkblue",lty=2,lwd=2)

  Y <- c(dy.dx$y.fit[1],dy.dx$y.fit[1]+cumsum(diff(dy.dx$y.fit)))
  #print(Y)
  lines(X,Y,col="grey",lty=2)
}
