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
}
\keyword{manip}


