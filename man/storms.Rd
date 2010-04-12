\name{Storms.ERA40}
\alias{Storms.ERA40}
\alias{Storms.ERA40x1}
\alias{Storms.NMC}
\alias{data(Storms.NMC)}
\alias{data(Storms.ERA40)}
\title{North Atlantic cyclone statistics.}
\description{Counts of low-pressure systems every 12hr over the North
  Atlantic region (00.00 and 12.00 UTC), derived from the NMC (now
  National Centers for Environmental Prediction, NCEP,
  \url{http://www.ncep.noaa.gov/}) or ECMWF re-analysis ERA40 6-hr or 24-hr slp
  fields ('x1' is once a day).}
\usage{data(Storms.ERA40)
       data(Storms.ERA40x1)
       data(Storms.NMC) }
\format{'Storms.ERA40' - a list of lon (longitude in degrees), lat (latitude in
  degrees), psl (central sea level pressure in hPa), yy (year), mm (month), dd (day),
  max.dpsl (max pressure gradient), tim (time index), 
  max.speed (estimated maximum geostropic speed), radius (radius of low
  pressure system in km), rad.max.dpsl, dx, dy, version, scale.factor.x, scale.factor.y.
        'storms' - a list of lon, lat, psl, yy, mm, dd. See \code{\link{CCI}}.}
\source{The Norwegian Meteorological Institute, Climatology division.}
\references{EMS04-A-00146, Benestad, R.E.; Sorteberg, A.; Chen, D. 'Storm statistics derived using a calculus-based cyclone identification method', \url{http://www.cosis.net/members/meetings/sessions/oral_programme.php?p_id=110&s_id=1845}, European Meteorological Society AC2, Nice, Sept 28, 2004.}
\examples{
data(Storms.NMC)
summary(Storms.NMC)
cyclstat(Storms.NMC)
}
\keyword{datasets}
