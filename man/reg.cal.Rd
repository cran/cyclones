\name{Regression-gradients}
\alias{coefFit}
\alias{coefDeriv}
\alias{coefInt}
\alias{geoGrad}
\alias{testReg.cal}
\alias{derivFFT}
\alias{integrFFT}
\alias{testFFTcalc}
\title{ REGression based CALculus for empirical data}
\description{ Use an n-order polynomial fit to model a data series. Ths polynomial
 then provides the basis for the calculus: derivations and integrations.
 One way of estimating slopes in the terrain or the geostrophic winds
 components from sea level pressure fields.

\code{coefFit} fits the coeficients of a power series \deqn{y= a_0 +
a_1 x + a_2 x^2 + a_3 x^3 + ...}{y= a0 + a1 x + a2 x^2 + a3 x^3 + ...}, where \deqn{x = -1,...,+1}{x = -1,...,+1}. Reference: R.E. Benestad (2003) What can present climate models tell us about climate change? Climatic Change Vol 59, 311-332

\code{coefDeriv} computes the first derivative of a power series \deqn{y= a_0 +
a_1 x + a_2 x^2 + a_3 x^3 + ...}{y= a0 + a1 x + a2 x^2 + a3 x^3 + ...} using the results from \code{coefFit}  so that \deqn{dy/dx= a_1  + 2 a_2 x
+ 3 a_3 x^2 + ...}{dy/dx= a1 + 2 a2 x + 3 a3 x^2 +
...}. (Reference:e.g. G. Stephenson (1961), Mathematical methods for
science students, Logman  Scientific \& Technical, p. 86).

\code{coefInt} integrates a series \deqn{y= a_0 + a_1 x + a_2 x^2
+ a_3 x^3 + ...}{y= a0 + a1 x + a2 x^2 + a3 x^3 + ...} using the
results from \code{coefFit} so that \deqn{\int ydx = a_0 x + 1/2 a_1
x^2 + 1/3 a_2 x^3 + 1/4 a_3 x^4 +  ... + const}{integral(y)=a0 x + 1/2 a1
x^2 + 1/3 a2 x^3 + 1/4 a3 x^4 +  ... + const}.

\code{geoGrad} computes gradients of the topography.The topolography data
file can be obtained from \url{http://ferret.pmel.noaa.gov/NVODS/servlets/dataset}.

\code{testReg.cal} is a test-function for \code{coefFit},
\code{coefDeriv} and \code{coefInt}.

\code{derivFFT} uses FFT to compute the first derivative. Anolagous to \code{\link{dX}}

\code{integrFFT}  uses FFT to integrate a series and is the (pseudo) inverse of
\code{derivFFT}. 

\code{testFFTcalc} is a test function for \code{derivFFT} and \code{integrFFT}.
} 
\usage{
coefFit(y,x=NULL,n=length(y),method="lm")
coefDeriv(y)
coefInt(y,c1=0)
geoGrad(maxhar = 35, fname = "data/etopo5_scandinavia.Rdata",
        x.rng = c(-10, 35), y.rng = c(50, 73), plot = FALSE)
testReg.cal(i.y=240,N=50)
derivFFT(y)
integrFFT(y)
testFFTcalc()
}
\arguments{
    \item{y}{A vector.}
    \item{x}{A vector of index values (e.g. time or longitude).}
    \item{n}{Number of harmonics to fit}
    \item{method}{Regression model for fitting the harmonics}
    \item{c1}{first coefficient: the constant term}
    \item{maxhar}{number of harmonix (see \code{\link{dX}})}
    \item{i.y}{}
    \item{N}{}
    \item{x.rng}{longitude range}
    \item{y.rng}{latitude range}
    \item{fname}{Name of topology file}
    \item{plot}{Plot the results}
  }
\value{}
\author{R.E. Benestad}
\examples{
\dontrun{
# Polynomial series calculus. (testReg.cal)
  library(clim.pact)
  load("data/etopo5_scandinavia.Rdata")
  y <- ROSE[i.y,]
  a <- coefFit(y,n=N)
  da <- coefDeriv(a)
  a.2 <- coefInt(da,c1=a$coefs[1])
  newFig()
  plot(ETOPO5X,y,type="s",lwd=3,xlab="Longitude (degE)",ylab="m.a.s.l.",
       main=paste("Transect: ",round(ETOPO5Y[i.y],1),"degE"))
  polygon(c(ETOPO5X,ETOPO5X[1320],ETOPO5X[1]),
          c(da$y.deriv/quantile(da$y.deriv,0.9)*quantile(y,0.7),0,0),col="blue")
  lines(ETOPO5X,y,type="s",lwd=4)
  lines(ETOPO5X,a$y.hat,col="red",lty=2,lwd=2)
  lines(ETOPO5X,a.2$y.int,col="steelblue",lty=1)


# FFT-based calculus  (testFFTcalc)

x <- seq(-3,3,length=100)
y <- 3*cos(5*x) + 0.3*sin(18*x) - 1.4*sin(3*x) + 0.15*sin(23*x)
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
}
\keyword{manip}


