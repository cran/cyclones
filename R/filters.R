# A number of different filters using different window
# shapes.
#
# R.E. Benestad, July, 2002, met.no.
#
# ref: Press et al. (1989), Numerical Recipes in Pascal, pp. 466
#library(ts)

# Moving-average (box-car) filter
ma.filt <- function(x,n) {
  y <- filter(x,rep(1,n)/n)
  y
}

# Gaussian filter with cut-off at 0.025 and 0.975 of the area.
gauss.filt <- function(x,n) {
  i <- seq(0,qnorm(0.975),length=n/2)
  win <- dnorm(c(sort(-i),i))
  win <- win/sum(win)
  y <- filter(x,win)
  y
}

# Binomial filter
binom.filt <- function(x,n) {
  win <- choose(n-1,0:(n-1))
  win <- win/max(win,na.rm=T)
  win[is.na(win)] <- 1
  win <- win/sum(win,na.rm=T)
  y <- filter(x,win)
  y
}

# Parzen filter (Press,et al. (1989))
parzen.filt  <-  function(x,n) {
  j <- 0:(n-1)
  win <- 1 - abs((j - 0.5*(n-1))/(0.5*(n+1)))
  win <- win/sum(win)
  y <- filter(x,win)
  y
}

# Hanning filter (Press,et al. (1989))
hanning.filt  <-  function(x,n) {
  j <- 0:(n-1)
  win <- 0.5*(1-cos(2*pi*j/(n-1)))
  win <- win/sum(win)
  y <- filter(x,win)
  y
}

# Welch filter (Press,et al. (1989))
welch.filt  <-  function(x,n) {
  j <- 0:(n-1)
  win <- 1 - ((j - 0.5*(n-1))/(0.5*(n+1)))^2
  win <- win/sum(win)
  y <- filter(x,win)
  y
}

