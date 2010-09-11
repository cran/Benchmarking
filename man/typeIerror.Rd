% $Id: typeIerror.Rd 72 2010-09-11 17:06:14Z Lars $
\name{typeIerror}
\alias{typeIerror}

\title{Probability of a type I error for a test in bootstraped DEA
models}

\description{Calculates the probability of a type I error for a test
in bootstraped DEA models; an after method for the method boot.sw98 in
the FEAR package}

\usage{typeIerror(shat,s) }

\arguments{

\item{shat}{The value of the statistic for which the propability of a
type I error is to be calculated}

\item{s}{Vector with calculated values of the statistic for each of
the \code{NREP} bootstraps; NREP is from \code{boot.sw98} }

}

\details{Needs bootstraped values of the test statistic}

\value{Returns the probability of a type I error}



\author{Lars Otto \email{lo@foi.dk}}


\examples{

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
