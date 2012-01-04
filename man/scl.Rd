\name{Benestad2005}
\alias{sunspots}
\alias{f10.7cm}
\alias{Benestad2005}
\alias{gcr}
\title{Analysis of solcar cycle length (SCL).}
\description{
The R-code used for making the analysis and figures in Benestad,
R.E. (2005) 'A review of the solar cycle length estimates' GRL 32
L15714, doi:10.1029/2005GL023621, August 13. Note. the 10.7 cm flux used
to be fetched from URL
ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/DAILYPLT.OBS
which no longer seems to work.
}

\usage{
f10.7cm(url="ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/Penticton_Adjusted/daily/DAILYPLT.ADJ",plot=TRUE,na.rm=TRUE)
sunspots(url="http://sidc.oma.be/DATA/monthssn.dat",plot=TRUE)
Benestad2005()
data(gcr)
}
\arguments{
    \item{url}{URL}
    \item{plot}{Option to produce graphics}
    \item{na.rm}{TRUE: Ignore NAs}
}
\value{A table or lists containing the relevant data.}
\author{R.E. Benestad}
\examples{
Benestad2005()
}
\keyword{manip}


