\name{within}
\alias{within}
\alias{within.old}
\alias{test.within}
\title{Points within an area.}
\description{
Finds points within a given closed 2-D boundary.
}

\usage{
within(x,y,X,Y,plot=FALSE,test=FALSE)
within.old(x,y,X,Y,test=FALSE)
test.within(a=3,b=5)
}
\arguments{
    \item{x}{vector with x-cordinates to test if they lie inside given boundary.}
    \item{y}{vector with y-cordinates to test...}
    \item{X}{vector of x-coordinates of the boundary.}
    \item{Y}{vector of y-coordinates of the boundary.}
    \item{plot}{TRUE if plotting results.}
    \item{test}{TRUE if conducting test.}
    \item{a}{x-axis of ellipse.}
    \item{b}{y-axis of ellipse.}    
  }
\value{Bolean vector of same length as \code{x}.}
\author{R.E. Benestad}
\examples{
x <- rnorm(100)
y <- rnorm(100)
X <- cos(2*pi*seq(0,1,length=360))
Y <- cos(2*pi*seq(0,1,length=360))
print(sum(within(x,y,X,Y)))
}
\keyword{manip}


