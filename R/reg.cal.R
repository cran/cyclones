# REGression based CALculus for empirical data
#----------------------------------------------
# R.E. Benestad.
# Use an n-order polynomial fit to model a data series. This polynomial
# then provides the basis for the calculus: derivations and integrations.
# One way of estimating slopes in the terrain or the geostrophic winds
# components from sea level pressure fields.
#
# Addition: using FT.
#
# Ref. R.E. Benestad,
# "What can present climate models tell us about climate change?"
# Climatic Change, accepted 2002.
#
# REB: Oslo, November 5 2002.
# also see dx.R &  dy.R

coefFit <- function(y,x=NULL,n=length(y),method="lm") {
  if (is.null(x))  x <- seq(-1,1,length=length(y))
  exprn <- paste(method,"(y ~ 1",sep="")
  for (i in 1:n) {
    exprn <- paste(exprn," + I(x^",i,")",sep="")
  }
  exprn <- paste(exprn,")",sep="")
#  print(exprn)
  expr <- parse(text=exprn)
  n.ordr.lm <- eval(expr)
  stat <- summary(n.ordr.lm)
#  print(stat)
  y.hat <- predict(n.ordr.lm)
  coef.fit <- list(coefs=as.real(n.ordr.lm$coefficients),
                   y.hat=y.hat,y=y,x=x,model=n.ordr.lm,n=n)
  class(coef.fit) <- "coef.fit"
  invisible(coef.fit)
}

coefDeriv <- function(y) {
  if (class(y)!="coef.fit") {
    stop("This function requires a 'coef.fit' object")
  }
  coef.deriv <- rep(0,(y$n+1))
  for (i in 1:(y$n+1)) {
    coef.deriv[i]  <-  (i-1)*y$coefs[i]
  }
  y.deriv <- y$y*0
  for (i in 1:(y$n+1)) {
    if (is.finite(coef.deriv[i])) y.deriv <- y.deriv + coef.deriv[i]*y$x^(i-1)
  }
  coef.der<- list(y.deriv=y.deriv,coef.deriv=coef.deriv,
                  coefs=y$coef,
                  y.hat=y$y.hat,y=y$y,x=y$x,model=y$model,n=y$n)
  class(coef.der) <- c("coef.fit","coef.deriv")
  invisible(coef.der)
}

coefInt <- function(y,c1=0) {
  if (class(y)[2]!="coef.deriv") {
    stop("This function requires a 'coef.fit' object")
  }
  coef.int <- rep(0,(y$n+1))
  coef.int[1] <- c1
  for (i in 2:(y$n+1)) {
    coef.int[i]  <-  y$coef.deriv[i]/(i-1)
  }
  y.int <- y$y*0
  for (i in 1:(y$n+1)) {
    if (is.finite(coef.int[i])) y.int <- y.int + coef.int[i]*y$x^(i-1)
  }
  coef.int<- list(y.int=y.int,coef.int=coef.int,
                  coefs=y$coef,
                  y.hat=y$y.hat,y=y$y,x=y$x,model=y$model,n=y$n)
  class(coef.int) <- c("coef.fit","coef.int")
  invisible(coef.int)
}

geoGrad <- function(maxhar=35,fname="data/etopo5_scandinavia.Rdata",x.rng=c(-10,35),y.rng=c(50,73),plot=FALSE) {
  library(clim.pact)
  if (file.exists(fname)) load(fname) else {
     print("Try:")
     print("URL: http://www.unidata.ucar.edu/cgi-bin/dods/datasets/datasets.cgi?keyword=etopo5&xmlfilename=datasets.xml")
     print("URL: http://ferret.pmel.noaa.gov/cgi-bin/dods/nph-dods/data/PMEL/etopo5.nc.html")
  }
  print("Do not reduce the region")

  if (!is.null(x.rng)) {   x.keep <- (ETOPO5X >= x.rng[1]) & (ETOPO5X <= x.rng[2]) }
  if (!is.null(y.rng)) {   y.keep <- (ETOPO5Y >= y.rng[1]) & (ETOPO5Y <= y.rng[2]) }
  ROSE <- ROSE[y.keep,x.keep]
  ETOPO5X <-ETOPO5X[x.keep]
  ETOPO5Y <-ETOPO5Y[y.keep]
  nxy <- dim(ROSE)
  grad.ew <- ROSE*0
  grad.ns <- ROSE*0

  if (plot) {
    newFig()
    image(ETOPO5X,ETOPO5Y,t(ROSE))
    addland()
    grid()
    dev.copy2eps(file="geoGrad_z.eps")
  }

  print("GRADIENT E-W")
  grad.ew <- dX(ETOPO5X,ETOPO5Y,t(ROSE),maxhar=maxhar,plot=plot)
  if (plot) {
    newFig()
    image(ETOPO5X,ETOPO5Y,grad.ew$dZ)
    addland()
    grid()
    dev.copy2eps(file="geoGrad_dz.dx.eps")
  }

  print("GRADIENT N-S")
  grad.ns <- dY(ETOPO5X,ETOPO5Y,t(ROSE),maxhar=maxhar,plot=plot)
  
  if (plot) {
    newFig()
    image(ETOPO5X,ETOPO5Y,grad.ns$dZ)
    addland()
    grid()
    dev.copy2eps(file="geoGrad_dz.dy.eps")
  }

  slope <- sqrt(grad.ew$dZ^2 + grad.ns$dZ^2)
  if (plot) {
    newFig()
    image(ETOPO5X,ETOPO5Y,slope)
    addland()
    grid()
    dev.copy2eps(file="geoGrad_slope.eps")
  }
  direction <- atan2(grad.ns$dZ,grad.ew$dZ)
  if (plot) {
    newFig()
    image(ETOPO5X,ETOPO5Y,direction)
    addland()
    grid()
    dev.copy2eps(file="geoGrad_direction.eps")
  }
  save(file="data/geo.grad2.Rdata",grad.ew,grad.ns,
       slope,direction,ETOPO5X,ETOPO5Y)
}

testReg.cal <- function(i.y=240,N=50) {
  load("data/etopo5_scandinavia.Rdata")
  y <- ROSE[i.y,]
  a <- coefFit(y,n=N)
  da <- coefDeriv(a)
  a.2 <- coefInt(da,c1=a$coefs[1])
#  print(a$coefs)
#  print(c(0,a$coefs[2],2*a$coefs[3],3*a$coefs[4]))
#  print(da$coef.deriv)
#  print(c(a$coefs[1],da$coef.deriv[2],da$coef.deriv[3]/2,da$coef.deriv[4]/3))
#  print(a.2$coef.int)
  newFig()
  plot(ETOPO5X,y,type="s",lwd=3,xlab="Longitude (degE)",ylab="m.a.s.l.",
       main=paste("Transect: ",round(ETOPO5Y[i.y],1),"degE"))
  polygon(c(ETOPO5X,ETOPO5X[1320],ETOPO5X[1]),
          c(da$y.deriv/quantile(da$y.deriv,0.9)*quantile(y,0.7),0,0),col="blue")
  lines(ETOPO5X,y,type="s",lwd=4)
  lines(ETOPO5X,a$y.hat,col="red",lty=2,lwd=2)
  lines(ETOPO5X,a.2$y.int,col="steelblue",lty=1)
  grid()
}

derivFFT <- function(y) {
# P(0) = 1/N^2 |C_0|^2
# P(f_k) = 1/N^2 (|C_k|^2 + |C_{N-k}|^2)
# P(f_N) = 1/N^2 |C_{N/2}|^2
# Y  = a0 + a1 cos(t w1) + b1 sin(t w1)  +  a2 cos(t w2) + b2 sin(t w2) + ...
# Y' =  -a1/w1 sin(t w1) + b1/w1 cos(t w1) - a2/w2 sin(t w2) + b2/w2 cos(t w2) + ...
# Im(Y) -> sin; Re(Y) -> cos
  Y <- fft(y,inverse=FALSE)
  n <- length(Re(Y))
  Y.im <- -1*Re(Y)
  Y.re <- Im(Y)
  Y[1] <- 0
  print(c(n,n/2))
  if (n/2 - trunc(n/2)==0) {
    for (i in 2:(n/2-1)) {
      Y.im[i] <- Y.im[i]*(i-1)
      Y.re[i] <- Y.re[i]*(i-1)
      Y.im[n-i+2] <- Y.im[n-i+2]*(i-1)
      Y.re[n-i+2] <- Y.re[n-i+2]*(i-1)
    }
    Y.im[n/2] <- Y.im[n/2]*n/2
    Y.re[n/2] <- Y.re[n/2]*n/2
  } else {
    for (i in 2:(trunc(n/2)+1)) {
      Y.im[i] <- Y.im[i]*(i-1)
      Y.re[i] <- Y.re[i]*(i-1)
      Y.im[n-i+2] <- Y.im[n-i+2]*(i-1)
      Y.re[n-i+2] <- Y.re[n-i+2]*(i-1)
    }
  }
#  lines(seq(-3,3,length=100),Im(Y)/sd(Im(Y)),col="darkblue",type="s")
#  lines(seq(-3,3,length=100),Re(Y)/sd(Re(Y)),col="darkred",type="s")
  Y <- complex(real = Y.re, imaginary = Y.im)
#  lines(seq(-3,3,length=100),Y.re/sd(Y.re),col="steelblue",type="s",lwd=2,lty=2)
#  dydx <- Re(fft(Y,inverse=TRUE)/length(y))
  dydx <- fft(Y,inverse=TRUE)/length(y)
  invisible(dydx)
}

integrFFT <- function(y) {
# P(0) = 1/N^2 |C_0|^2
# P(f_k) = 1/N^2 (|C_k|^2 + |C_{N-k}|^2)
# P(f_N) = 1/N^2 |C_{N/2}|^2
# Y  = a0 + a1 cos(t w1) + b1 sin(t w1)  +  a2 cos(t w2) + b2 sin(t w2) + ...
# Y' =  -a1/w1 sin(t w1) + b1/w1 cos(t w1) - a2/w2 sin(t w2) + b2/w2 cos(t w2) + ...
# Im(Y) -> sin; Re(Y) -> cos
  Y <- fft(y,inverse=FALSE)
  n <- length(Re(Y))
  Y.im <- Re(Y)
  Y.re <- -1*Im(Y)
  Y.im[1] <- 0
  Y.re[1] <- 0
  if (n/2 - trunc(n/2)==0) {
    for (i in 2:(n/2-1)) {
      Y.im[i] <- Y.im[i]/(i-1)
      Y.re[i] <- Y.re[i]/(i-1)
      Y.im[n-i+2] <- Y.im[n-i+2]/(i-1)
      Y.re[n-i+2] <- Y.re[n-i+2]/(i-1)
    }
    Y.im[n/2] <- Y.im[n/2]/(n/2)
    Y.re[n/2] <- Y.re[n/2]/(n/2)
  } else {
    for (i in 2:(trunc(n/2)+1)) {
      Y.im[i] <- Y.im[i]/(i-1)
      Y.re[i] <- Y.re[i]/(i-1)
      Y.im[n-i+2] <- Y.im[n-i+2]/(i-1)
      Y.re[n-i+2] <- Y.re[n-i+2]/(i-1)
    }
  }
#  lines(seq(-3,3,length=100),Im(Y)/sd(Im(Y)),col="darkblue",type="s")
#  lines(seq(-3,3,length=100),Re(Y)/sd(Re(Y)),col="darkred",type="s")
  Y <- complex(real = Y.re, imaginary = Y.im)
#  lines(seq(-3,3,length=100),Y.re/sd(Y.re),col="steelgrey",type="s",lwd=2,lty=2)  
#  y.int <- Re(fft(Y,inverse=TRUE)/length(y))
  y.int <- fft(Y,inverse=TRUE)/length(y)
  invisible(y.int)
}  

testFFTcalc <- function() {
  x <- seq(-3,3,length=100)
#  y <- 1.4*x + 0.7*x^2 - 0.3*x^3 + 0.1*x^4 - 0.05*x^5 + 0.01*x^6
  y <- 3*cos(5*x) + 0.3*sin(18*x) - 1.4*sin(3*x) + 0.15*sin(23*x)
  print(summary(y))
  plot(x,y,type="l",lwd=3)
  grid()
  dydx <- deriv.fft(y)
  y2 <- integr.fft(dydx)
  lines(x,y2,col="red",lty=2,lwd=2)
  grid()
  newFig()
  plot(x,Im(dydx),type="l",lwd=3)
  grid()
  lines(x,rep(0,length(x)),col="grey")  
}
