\name{cyclstat}
\alias{cyclstat}
\title{Cyclone statistics}
\description{Produce storm statistics/analysis. This function was used to produce the figures in the paper Benestad, Sorteberg and Chen, 'The use of a Calculus-based Cyclone Identification method for generating storm statistics' submitted to Climate Dynamics 2004.}
\usage{
cyclstat(fname=NULL,psl0=1000,topo="data/etopo60.Rdata",cyclone=TRUE,
         x.rng=c(5,35),y.rng=c(55,72),cmp=FALSE,mon=c(12,1,2))
}
\arguments{
    \item{fname}{Filename of cyclone stats data.}
    \item{psl0}{Treshold value}
    \item{topo}{Name of topography file.}
    \item{cyclone}{Treshold: TRUE -> less than, FALSE -> greater than}
    \item{x.rng}{Interval of longitudes of interest (deg E)}
    \item{y.rng}{Interval of latitudes of interest (deg N)}
    \item{cmp}{Flag set for comparing two storm data sets: one given a fname and the default in cyclones ('data(storms)')}
    \item{mon}{Months to plot a histogram for.}
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


