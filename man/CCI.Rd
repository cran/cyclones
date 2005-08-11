\name{CCI}
\alias{CCI}
\alias{stopCCI}
\title{Calculus Cyclone identification.}
\description{
Identifies cyclones (low pressure systems) in a gridded data set using a Calculus Cyclone Identification (CCI) method (EMS04-A-00146, Benestad, R.E.; Sorteberg, A.; Chen, D. 'Storm statistics derived using a calculus-based cyclone identification method', \url{http://www.cosis.net/members/meetings/sessions/oral_programme.php?p_id=110&s_id=1845}, European Meteorological Society AC2, Nice, Sept 28, 2004). Storms are identified with longitude, latitude, and date. Also returned are estimates of local minimum pressure, max pressure gradient near storm, max geostrophic windspeed near storm, and radius of the storm. The storm location is by means of finding where first derivatives of north--south and east--west gradients both are zero. The storm radius is estimated from the points of inflexion along the latitude and longitude lines running trough the centre of the storm. 

This code is the basis for the analysis carried out by Benestad et al. 'The use of a Calculus-based Cyclone Identification method for generating storm statistics' submitted to Climate Dynamics (2004).

\code{stopCCI} signals the CCI process to stop in a safe manner by removing the file '.CCI.run' in the run directory. 
} 
\usage{
CCI(maxhar=25,lplot=TRUE,nsim=10,fname="data/cyclones.Rdata",
    fielddata="data/nmc_slp.nc",vname="slp",cyclones=TRUE,force365.25=FALSE,
    x.rng=c(-80,40),y.rng=c(20,75),tslice=3652,rad=5,dx=1,dy=1,
    times=NULL,label=NULL,rho=1.293,nc.read.method="retrieve.nc")
stopCCI()
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
    \item{dx}{Resolution in longitude (in degrees).}
    \item{dy}{Resolution in latitudes (in degrees).}
    \item{times}{Specify time slices for reading chunks of data.}
    \item{label}{Label for ID-purposes.}
    \item{rho}{Density of air for calculation of geostrophic wind.}
    \item{nc.read.method}{Method to read netCDF or other gridded data. Can be set to a custom made routine.THe output should be a 'field'-type object.}
  }
\value{
A list: \code{list(lon=lon, lat, tim, psl, yy, mm, dd, i, label, max.dpsl, max.speed, radius, rad.max.dpsl, dx, dy)}. The subobjects \code{lon} (longitude: units degrees), \code{lat} (latitude: units: degrees), \code{psl} (local minimum pressure: units hPa), \code{max.dpsl} (pressure gradient: units Pa/m), \code{max.speed} (windspeed: units m/s), and \code{radius} (radius of the storm: units km) are \code{[1:nt,1:i.sim]}-matrices. 
}
\author{R.E. Benestad}}
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




