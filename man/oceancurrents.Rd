\name{FrontStats}
\alias{FrontStats}
\alias{Gulfstream}
\alias{Kuroshio}
\alias{CurrentClimatology}
\alias{CurrentHistory}
\alias{CurrentSpectrum}
\alias{Current2monthly}
\alias{CurrentLatitude}
\alias{CurrentValidation}
\alias{data(gulfstream)}
\alias{data(kuroshio)}
\title{Ocean currents and sea surface temperature fronts}
\description{
The function \code{FrontStats} uses \code{dY} to identify the
latitudinal position of the sharpest stemperature (sea sufrace
temperature=SST) gradient (dT/dy) along a series
of north-south transects. \code{Gulfstream} and \code{Kuroshio} use
\code{FrontStats} to generate a set of coefficients describing the shape
and position of the Gulfstream and the Kuroshio ocean currents, based on
the position of the sharpest latitudinal SST gradient (dT/dy). The position and
the shape of the ocean current - the sharpest latitudinal SST gradient -
is described approximately by equation

\deqn{\hat{y} = \beta_0 +  \beta_1 x + \beta_2 x^2 +  \beta_3 x^3 +
  \beta_4 x^4 + \beta_5 x^5 + \beta_6 x^6}{y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6)}.

The best-fit solution is derived using multiple linear regression (\code{lm}).

A demonstration/validation of these results is available as an animation
at YouTube: \url{http://www.youtube.com/watch?v=SQJBa7ynUpw} and \url{http://www.youtube.com/watch?v=FgDgKhgJ5WA}

}

\usage{
FrontStats(sst,plot=FALSE,plot.int=100,maxhar=6,neg.gradient=TRUE)
Gulfstream(fname="sst.wkmean",plot=FALSE,maxhar=12,plot.int = 100)
Kuroshio(fname="sst.wkmean",plot=FALSE,maxhar=12,plot.int = 100)
CurrentClimatology(current="gulfstream",parameter=1) 
CurrentHistory(current="gulfstream",parameter=1)
CurrentSpectrum(current="gulfstream",parameter=1,frequency=FALSE)
Current2monthly(current="gulfstream",parameter=1)
CurrentLatitude(current="gulfstream",lon=NULL,nx=50,plot=TRUE)
CurrentValidation(current="gulfstream",fname="sst.wkmean",plot=TRUE,las=1)
data(gulfstream)
data(kuroshio)
}
\arguments{
    \item{sst}{A field object (here the sea surface temperature; SST)}
    \item{plot}{Flag for plotting}
    \item{plot.int}{Plotting interval}
    \item{maxhar}{Number of harmonics to include. See \code{\link{dX}}.}
    \item{neg.gradient}{Look for negative gradients (decreasing SST with latitude)}
    \item{fname}{Name pattern of netCDF files with SST.}
    \item{lon}{a range of longitude coordinates defining reion of interest}
    \item{nx}{length of }
    \item{las}{as \code{las} in \code{\link{par}}}
    \item{current}{Name of current system ("gulfstream"/"kuroshio")}
    \item{parameter}{parameter/coefficient of the polynomial fit to SST front}
    \item{frequency}{plot spectrum against frequency or periodicity.}
      
}
\value{
A matrix containing the 9 regression coefficients, the R-squared value
and the 95-percentile of dT/dy.  
}
\author{R.E. Benestad}
\examples{
\dontrun{
sst <- retrieve.nc("sst.wkmean.1981-1989.nc",x.rng=c(-80,-30),y.rng=c(30,50),v.nam="sst") 
coefs <- FrontStats(sst,plot=plot,maxhar=16,plot.int=1)

Gulfstream() -> gulfstream; save(file="gulfstream.rda",gulfstream)

Kuroshio() -> kuroshio; save(file="kuroshio.Rdata",kuroshio)

# This call only works completely for Linux platform with ImageMagic and
# gifmerge installed. Type 'which convert' and 'which gifmerge' to see
# if they are installed. These can be obtained from:
# http://www.imagemagick.org/
# URL: http://the-labs.com/GIFMerge/
# 
CurrentValidation()
}
}
\keyword{manip}


