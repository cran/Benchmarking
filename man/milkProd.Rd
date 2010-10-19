% $Id: milkProd.Rd 80 2010-10-20 08:58:43Z Lars $
\name{milkProd}
\alias{milkProd}
\docType{data}

\title{ Data: Milk producers }

\description{
Data colected from Danish  milk producers.
}

\usage{data(milkProd)}

\format{
  A data frame with 108 observations on the following 5 variables.
  \describe{
    \item{\code{farmNo}}{farm number}
    \item{\code{milk}}{Output of milk, kg}
    \item{\code{energy}}{Energy expenses}
    \item{\code{vet}}{Veterinary expenses}
    \item{\code{cows}}{Number of cows}
  }
}

% \details{
% }

\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header = TRUE, sep = ";")} such that this file
is a semicolon separated file and not a comma separated file.  }

\source{
Accounting and buissiness check data
}

\examples{
data(milkProd)
pairs(milkProd,gap=0)
}
\keyword{datasets}
