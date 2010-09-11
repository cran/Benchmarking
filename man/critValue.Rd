% $Id: critValue.Rd 72 2010-09-11 17:06:14Z Lars $
\name{critValue}
\alias{critValue}

\title{Critical values from bootstraped DEA models}

\description{Calculates critical value for test using bootstrap output
in DEA models; an after method for the method boot.sw98 in the FEAR
package} 

\usage{ critValue(s,alfa=0.05) }

\arguments{

\item{s}{Vector with calculated values of the statistic for each of
the \code{NREP} bootstraps; NREP is from \code{boot.sw98} }

\item{alfa}{The size of the test}

}

\details{Needs bootstraped values of the test statistic}

\value{Returns the critical value}


\author{Lars Otto \email{lo@foi.dk}}


\examples{


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
