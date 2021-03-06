% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploreRF.R
\name{exploreRF}
\alias{exploreRF}
\alias{plot.RFw}
\title{Plot ....}
\usage{
exploreRF(rf, somC, data, k = 10, ...)
}
\arguments{
\item{rf}{Random forest model resulted from \code{\link{wrapRF}}}

\item{somC}{object resulted from \code{\link{wrapSOM}} function}

\item{data}{numerical data.frame}

\item{k}{Number of variables to plot}
}
\value{
The function returns ...
}
\description{
Plot results from \code{RFwrap} ...
}
\details{
...
}
\examples{
data(sim1)
rf=wrapRF(nema_PCRBS,somC)
exploreRF(rf, somC, nema_PCRBS)
}
\seealso{
...
}
\author{
Danilo Candido Vieira
}
