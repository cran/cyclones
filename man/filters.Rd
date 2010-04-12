\name{filters}
\alias{ma.filt}
\alias{gauss.filt}
\alias{binom.filt}
\alias{parzen.filt}
\alias{hanning.filt}
\alias{welch.filt}
\title{Various filters.}
\description{
A bundle of window types. See e.g. Press et al. (1989), Numerical Recipes in Pascal, pp. 466. i Moving-average (box-car) filter, ii Gaussian filter with cut-off at 0.025 and 0.975 of the area, iii Binomial filter, iv Parzen filter (Press,et al. (1989)), v Hanning filter (Press,et al. (1989)), vi Welch filter (Press,et al. (1989)).
}

\usage{
ma.filt(x,n)
gauss.filt(x,n)
binom.filt(x,n)
parzen.filt(x,n)
hanning.filt(x,n)
welch.filt(x,n)

}
\arguments{
    \item{x}{series to filter}
    \item{n}{window width}
  }
\value{A vector}
\author{R.E. Benestad}
\examples{
N <- 1000
t <- seq(1,N,by=1) * pi/180
x <- sin(2*t) + 0.3*rnorm(N)
plot(x)
lines(ma.filt(x,15),col="grey40",lwd=4)
lines(gauss.filt(x,15),col="red",lwd=4,lty=2)
lines(binom.filt(x,15),col="blue",lwd=2)
lines(parzen.filt(x,15),col="darkgreen",lwd=2,lty=2)
lines(hanning.filt(x,15),col="darkred")
lines(welch.filt(x,15),col="darkblue",lty=2)
legend(0,-1.2,c("ma.filt(x,15)","gauss.filt(x,15)","binom.filt(x,15)",
       "parzen.filt(x,15)","hanning.filt(x,15)","welch.filt(x,15)"),cex=0.7,
       lty=c(1,2,1,2,1,2),lwd=c(4,4,2,2,1,1),bg="grey95",
       col=c("grey40","red","blue","darkgreen","darkred","darkblue"))
}
\keyword{manip}


