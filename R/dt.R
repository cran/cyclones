# Calculates the t-derivatives for a time series
#
# R.E. Benestad, 27.04.2004.
#
# also see reg.cal.R, dx, dy

dT <- function(y,maxhar=NULL) {
  Y <- y
  nt <- length(y)
  if (is.null(maxhar)) maxhar <- nt
  maxhar <- min(nt,maxhar)
  W <- 2*pi/nt
  a <- rep(0,nt)
  b <- rep(0,nt)
  c <- rep(0,nt)
  dy <- rep(0,nt)
  y.fit <- rep(0,nt)
    for (iw in 1:maxhar) {
     good <- is.finite(y)
     if (sum(good)>10) {
      wt <- iw*W*seq(1,nt,by=1)
      x1 <- cos(wt); x2 <- sin(wt)
      harmfit <- data.frame(y=y, x1=x1, x2=x2)
      harmonic <- lm(y ~ x1 + x2,data=harmfit)
      y[good] <- harmonic$residual
      c[iw] <- harmonic$coefficients[1]
      a[iw] <- harmonic$coefficients[2]
      b[iw] <- harmonic$coefficients[3]
      if (!is.finite(c[iw])) c[iw] <- 0
      if (!is.finite(b[iw])) b[iw] <- 0
      if (!is.finite(a[iw])) a[iw] <- 0

      y.fit <- y.fit + a[iw]*cos(wt) + b[iw]*sin(wt) + c[iw]
      dy <- dy +iw*W*( -a[iw]*sin(wt) + b[iw]*cos(wt) )
    }
  }

  results <- list(y=Y,a=a,b=b,c=c,dy=dy,y.fit=y.fit)
  class(results) <- "dydt"
  attr(results,"long_name") <- "t-derivative"
  attr(results,"descr") <- "dT.R"
  invisible(results)
}
