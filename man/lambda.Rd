% $Id: lambda.Rd 238 2021-07-11 16:12:17Z X022278 $
\name{lambda}
\alias{lambda}
\alias{lambda.print}

\title{Lambdas or the weight of the peers }


\description{
The lambdas, i.e. the weight of the peers, for each firm.
}


\usage{
lambda(object, KEEPREF = FALSE)
lambda.print(x, KEEPREF = FALSE, ...)

}


\arguments{

  \item{object,x}{A Farrell object as returned from \code{\link{dea}} et
      al.}

  \item{KEEPREF}{if \code{TRUE} then all firms are kept as reference
      firms even though they have all zero weights (lambda); might
      come handy if one needs to calculate \code{X} x \code{lambda}
      such that the firms in \code{X} and \code{lambda} agree. If
      \code{FALSE}, the default, then only weight for the peers are in
      the matrix \code{lambda}. }

  \item{\dots}{Optional parameters for the print method.}


}


\details{ Only returns the lambdas for firms that appear as a
peer, i.e. only lambdas for firms where at least one element of the
lambda is positive.  }


\value{ The return is a matrix with the firms as rows and the peers as
columns.  }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}




\seealso{ \code{\link{dea}} }


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
