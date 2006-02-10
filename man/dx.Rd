\name{dX}
\alias{dX}
\alias{dY}
\alias{dT}
\alias{test.dX}
\alias{test.dY}
\alias{test.dT}
\title{First derivatives from harmonic fits.}
\description{

\code{dX} calculates the x-derivatives for gridded data in
longitude-latitude coordinates. Is based on the equations by Stephenson (1961) p. 162
\deqn{\frac{dA}{dx} = \frac{1}{r \cos(\phi)} * \frac{d}
{d\theta}}{dA/dx = 1/(r cos(PHI)) * d/dTHETA}
where PHI is the latitude in radians and THETA the longitude.

\code{dY} calculates the y-derivatives for gridded data in longitude-latitude coordinates. 
\deqn{\frac{dA}{dy} = \frac{1}{r} \frac{d}{d\phi}}{dA/dy= 1/r d/dPHI}.

\code{dT} calculates the e.g. time derivatives of a time series.

All functions solve for the gradients by fitting coeffiecients to a Fourier
approximation \deqn{p(x) = a_0 + \sum_i[a_i \cos(w_i x) + b_i \sin(w_i
x)]}{p(x) = a0 + sum(a[i] cos(w[i]*x) + b[i]sin(w[i]*x))} 
and dp(x)/dx is \deqn{\frac{dp(x)}{dx} = \sum[w_i(-a_i \sin(w_i x) +
b_i \cos(w_i x))]}{dp(x)/dx = sum(w[i](-a[i] sin(w[i]*x) + b[i]cos(w[i]*x)))}

Also see \code{\link{derivFFT}}.

Version 1.1-6:

(i) Correction to the wave number in \code{dX} and \code{dY} -
previously wave number was scaled as 1/[number of points] and gave incorrect magnitudes for the derivatives, but now the wave number is 1/distance.

(ii) Included a simple confidence check (chk.conf) to set coefficients
to zero when +- 1 sd includes zero. This check can be disabled. Note, it
is recommended to use the minimum value for maxhar.

(iii) Included new test routines to test the magnitudes of the gradients.
}

\usage{
dX(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE,plot=FALSE,chk.conf=1)
dY(lon,lat,Z,r=6.378e06,maxhar=NULL,mask.bad=TRUE,plot=FALSE,chk.conf=1)
dT(y,maxhar=NULL,plot=FALSE,chk.conf=1)
test.dX(field=NULL,maxhar=7)
test.dY(field=NULL,maxhar=5)
test.dT(y=NULL,maxhar=15)
}
\arguments{
    \item{lon}{Longitudes}
    \item{lat}{Latitude}
    \item{Z}{Field values of which the gradients are sought}
    \item{r}{mean radius of Earth}
    \item{maxhar}{Number of harmonics to include. If NULL, \code{maxhar} is set to the length of the series.}
    \item{mask.bad}{}
    \item{y}{a time series (vector).}
    \item{plot}{plot data and fit}
    \item{field}{Field to which the test is applied (default: DNMI.slp)}
    \item{chk.conf}{If not NULL, test whether conf. int. (+- chk.conf * sd) includes zero - if
                    so, set coefficient to zero.}
  }
\value{A list: list(Z=Z,a=a,b=b,c=c,dZ=dZ,Z.fit=Z.fit,lon=lon,lat=lat)
of class "map"}
\author{R.E. Benestad}
\examples{
# Construct a test series:
t <- seq(-3,3,by=0.1)
nt <- length(t)
y <- 2*rnorm(nt) + 10*sin(t/4) - 7*sin(t/3) + 4*cos(t) - sin(2*t) + 2*cos(3*t) + 0.5*sin(4*t)
plot(y)

dydt <- dT(y)
lines(dydt$y.fit,col="red")
dydt.2 <- dT(y,maxhar=3)
lines(dydt.2$y.fit,col="blue",lty=2)

# First derivative
plot(dydt$dy,type="l")
lines(dydt.2$dy,col="blue",lty=2)

# Second derivative
dy2dt2<-dT(dydt$dy)
dy2dt2.2<-dT(dydt.2$dy)
plot(dy2dt2$dy,type="l")
lines(dy2dt2.2$dy,col="blue",lty=2)
}
\keyword{manip}


