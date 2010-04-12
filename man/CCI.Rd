\name{CCI}
\alias{CCI}
\alias{stopCCI}
\alias{gradient.wind}
\title{Calculus Cyclone identification.}
\description{
Identifies cyclones (low pressure systems) in a gridded data set using a Calculus Cyclone Identification (CCI) method (EMS04-A-00146, Benestad, R.E.; Sorteberg, A.; Chen, D. 'Storm statistics derived using a calculus-based cyclone identification method', \url{http://www.cosis.net/members/meetings/sessions/oral_programme.php?p_id=110&s_id=1845}, European Meteorological Society AC2, Nice, Sept 28, 2004). Storms are identified with longitude, latitude, and date. Also returned are estimates of local minimum pressure, max pressure gradient near storm, max geostrophic windspeed near storm, and radius of the storm. The storm location is by means of finding where first derivatives of north--south and east--west gradients both are zero. The storm radius is estimated from the points of inflexion along the latitude and longitude lines running trough the centre of the storm. 

This code is the basis for the analysis carried out by Benestad et al. 'The use of a Calculus-based Cyclone Identification method for generating storm statistics' submitted to Climate Dynamics (2004).

If a north-south or east-west psl (sea level pressure) profile can be approximated as \deqn{p(\theta) = p_0 + \sum_{i=1}^{N_\theta} [ a_\theta(i) \cos(\omega_\theta(i) \theta) + b_\theta(i) \sin(\omega_\theta(i) \theta) ]}{p = p0 + sum [ a cos(w t) + b sin(w t) ]} 

Then the psl gradient can be estimated as: \deqn{\frac{\partial \hat{p}(\theta)}{\partial \theta} = \sum_{i=1}^{N_\theta} \omega_\theta(i)[ -\hat{a}_\theta(i) \sin(\omega_\theta(i) \theta) + \hat{b}_\theta(i) \cos(\omega_\theta(i) \theta)]}{dp/dt = sum (w[ -a sin(w t) + b cos(w w)])}

The gradient in x and y directions are found after the transform
\deqn{\frac{d\hat{p}(x)}{dx} = \frac{1}{a \cos(\phi)} \frac{d\hat{p}(\theta)}{d\theta}}{dp/dx= 1/[a cos(t)] dp(t)/dt} and \deqn{\frac{d\hat{p}(y)}{dy} = \frac{1}{a} \frac{d\hat{p}(\phi)}{d\phi}}{dp/dy = 1/a dp/dt} (Gill, 1982).



\code{stopCCI} signals the CCI process to stop in a safe manner by
removing the file '.CCI.run' in the run directory.

\code{gradient.wind} estimated the gradient wind of the data produced by
CCI after Fleagle and Businger (1980) p. 163. (eq 4.27). 

Reference: Benestad \& Chen (2006) 'The use of a Calculus-based Cyclone Identification method for generating storm statistics', Tellus A, in press.
           Benestad (2005) 
} 
\usage{
CCI(maxhar=14,lplot=TRUE,nsim=10,fname="data/cyclones.Rdata",
    fielddata="data/nmc_slp.nc",vname="slp",cyclones=TRUE,force365.25=FALSE,
    x.rng=c(-80,40),y.rng=c(5,75),tslice=3652,rad=5,dx=1,dy=1,
    times=NULL,label=NULL,rho=1.293,nc.read.method="retrieve.nc",
    graph.dir="CCI.graphs/",plot.interval=50,EPS=TRUE,verbose=TRUE,accuracy=NULL)
stopCCI()
gradient.wind(storms=NULL,icyclone=1)
}
\arguments{
    \item{maxhar}{number of harmonics used for fit to profile (Fourier truncation)}
    \item{lplot}{TRUE: produce plots.}
    \item{nsim}{Number of simultaneous cyclones identified and saved ordered according to depth/strength.}
    \item{fname}{Name of file containing the stats of cyclones (output).}
    \item{fielddata}{Name of file (netCDF) of gridded data (input)}
    \item{vname}{Variable name of the field in 'fielddata'.}
    \item{cyclones}{TRUE: identifies cyclones, FALSE anticyclones.}
    \item{force365.25}{TRUE - To ensure a 365.25-day year, FALSE evaluates to find the best match: 365.25-day or 360-day (model) year.}
    \item{x.rng}{Longitude region of interest.}
    \item{y.rng}{Latitude region of interest.}
    \item{tslice}{Number of maps read per chunk of data (time slice).}
    \item{rad}{Max radius of cyclone.}
    \item{dx}{Resolution in longitude (in degrees). NULL skips interpolation - uses original grid}
    \item{dy}{Resolution in latitudes (in degrees). NULL skips interpolation - uses original grid}
    \item{times}{Specify time slices for reading chunks of data.}
    \item{label}{Label for ID-purposes.}
    \item{rho}{Density of air for calculation of geostrophic wind.}
    \item{nc.read.method}{Method to read netCDF or other gridded
      data. Can be set to a custom made routine.THe output should be a
      'field'-type object.}
    \item{graph.dir}{Name of the directory in which to store the
      graphical output. Useful for checking the results, e.g. using
      \code{animate -delay 50 -pause 5 cci*.eps} under Linux.}
    \item{plot.interval}{Interval between the generation of each
      graphical output}
    \item{EPS}{TRUE for generation of encapsulated PostScript graphics, FALSE for PNG
      bitmaps. PNG bitmaps are useful for animatio (e.g. animate
      cci*.png; \url{http://www.imagemagick.org}), whereas postscript versions are good for hard copies.}
    \item{storms}{Data object returned by CCI with cyclone statistics.}
    \item{icyclone}{Decides which cyclone to do the gradient wind estimate for (eg 1 is for the deepest
      cyclone, 2 for the second deepest, etc).}
    \item{verbose}{Print out diagnostics.}
    \item{accuracy}{Not yet finished.}
  }
\value{
A list: \code{list(lon=lon, lat, tim, psl, yy, mm, dd, i, label,
  max.dpsl, max.speed, radius, rad.max.dpsl, dx, dy)}. The subobjects
\code{lon} (longitude: units degrees), \code{lat} (latitude: units:
degrees), \code{psl} (local minimum pressure: units hPa),
\code{max.dpsl} (pressure gradient: units Pa/m), \code{max.speed} 
(windspeed: units m/s if the units of SLP are in hPa), and \code{radius}
(radius of the storm: units km) are
\code{[1:nt,1:i.sim]}-matrices. \code{v.grad} is the estimated gradient wind.
}
\author{R.E. Benestad}
\examples{
\dontrun{
# Shell script for running in batch (background proccess)
#! /bin/bash
cat > paper20e.R << EOF
library(clim.pact)
library(cyclones)
source("cyclones/R/cyclones.R")

a<- Sys.info()
dir <- switch(substr(as.character(a[4]),1,9),
              "saragasso"="/data1/hirham/",
              "stratonim"="/data1/hirham/",
              "virvelvin"="/home/rasmusb/data/data1/")

filname.1 <- "/home/rasmusb/data/ERA40/era40_slp.nc"
fname.1 <- "data/cyclones_ERA40.Rdata"
vname.1 <- "msl"

print("===================<  ERA 15  >===================")
CCI(fname="data/cyclones_ERA40_r1.Rdata",dx=1,dy=1,fielddata=filname.1,
    vname=vname.1,label="ERA40: slp 1-degree res.",force365.25=TRUE,lplot=FALSE)
EOF

nice R CMD BATCH --no-save paper20e.R paper20e.out
}
}
\keyword{manip}




