\name{cyclstat}
\alias{cyclstat}
\alias{windstat}
\alias{discard.bad}
\title{Cyclone statistics}
\description{Produce storm statistics/analysis. This function was used
  to produce the figures in the paper Benestad, Sorteberg and Chen, 'The
  use of a Calculus-based Cyclone Identification method for generating
  storm statistics' accepted by Tellus A (2006). The function produces
  four plots: cyclstat.eps (a map showing frequency geographic
  distribution), cyclstat2.eps (time series showing evolution of
  number), cyclstat3.eps (annual cycle), and cyclstat4.eps (compare with
  Poisson distrubtion).
}
\usage{
cyclstat(fname=NULL,psl0=1000,topo="etopo60.Rdata",cyclone=TRUE,
         x.rng=c(5,35),y.rng=c(55,72),cmp=FALSE,mon=c(12,1,2),ERA40=TRUE)
windstat(fname=NULL,ws0=20,topo="etopo60.Rdata",cyclone=TRUE,
         x.rng=c(5,35),y.rng=c(55,72),cmp=FALSE,mon=c(12,1,2),ERA40=TRUE)
discard.bad(results)
}
\arguments{
    \item{fname}{Filename of cyclone stats data.}
    \item{psl0}{Treshold value: PSL lower than.}
    \item{ws0}{Treshold value: wind speed exceedance.}
    \item{topo}{Name of file containing topology; if doesn't exist, use
      data from the cyclone package.}
    \item{cyclone}{Treshold: TRUE -> less than, FALSE -> greater than}
    \item{x.rng}{Interval of longitudes of interest (deg E)}
    \item{y.rng}{Interval of latitudes of interest (deg N)}
    \item{cmp}{Flag set for comparing two storm data sets: one given a fname and the default in cyclones ('data(storms)')}
    \item{mon}{Months to plot a histogram for.}
    \item{ERA40}{True, use 6hr ERA40 as reference, otherwise use 12-hr
      NMC as reference.}
    \item{results}{CCI-object.}
  }
\value{}
\author{R.E. Benestad}
\examples{
\dontrun{
library(cyclones)
load("data/cyclonesERA40.Rdata")
results.era40 <- results

europe <- cyclstat(results.era40,cmp=TRUE,x.rng=c(-10,20),y.rng=c(40,55))
}
}
\keyword{manip}


