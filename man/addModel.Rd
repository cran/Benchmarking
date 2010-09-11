% $Id: addModel.Rd 72 2010-09-11 17:06:14Z Lars $
\name{add.dea}
\alias{add.dea}

\title{Additive DEA model}

\description{
Calculates efficiency as sum of input and output slacks
}


\usage{
add.dea(X, Y, RTS = "vrs", TRANSPOSE = FALSE, ...)
}


\arguments{

  \item{X}{input where efficiency is to be calculated, a K x m matrix
  of observations of K firms with m inputs (firm x input).  In case
  \code{TRABSPOSE=TRUE} the input matrix is transposed as input x
  firm. }

  \item{Y}{output where efficiency is to be calculated, a K x n matrix
  of observations of K firms with n outputs (firm x input).  In case
  \code{TRABSPOSE=TRUE} the output matrix is transposed as output x
  firm. }

  \item{RTS}{ A text string or a number to state the returns to scale
  assumption

\tabular{rcl}{
0 \tab "fdh" \tab Free disposability hull, no convexity assumption \cr
1 \tab "vrs" \tab Variable returns to scale, convex and free disposability \cr
2 \tab "drs" \tab Decreasing returns to scale, convex, down-scaling and free disposability \cr
3 \tab "crs" \tab constant returns to scale and free disposability \cr
4 \tab "irs" \tab increasing returns to scale, up-scaling, but not down-scaling \cr
5 \tab "add" \tab additivity,  scaling up and down, but only with integers
}

}

   \item{TRANSPOSE}{Input and output matrices are treated as goods times
      units as defaults.  If \code{TRANSPOSE=FALSE} then marices 
      are units times goods as in standard in R for statistical models}

   \item{\dots}{Optional parameters for the method \code{slack}.}

}


\value{

\item{sum}{Sum of all slacks for each unit}

\item{slack}{A non-NULL vector of logical variables, \code{TRUE} if
there is slack for the corresponding firm: also output from method
\code{slack}.}

\item{sx}{A matrix for input slacks for each firm, only calculated if
the option \code{SLACK} is \code{TRUE} or returned from the method
\code{slack}}

\item{sy}{A matrix for output slack, see \code{sx}}

}

\details{
 ~~ If necessary, more details than the __description__ above ~~
}

\source{
   corresponds to eqs. 4.34-4.38 in Cooper et al., 2007 
}

\references{
B&O
}

\examples{

x <- matrix(c(2,3,2,4,6,5,6,8),ncol=1)
y <- matrix(c(1,3,2,3,5,2,3,5),ncol=1)
dea.plot.frontier(x,y,txt=1:dim(x)[1])

sb <- add.dea(x,y,RTS="vrs")
data.frame("sx"=sb$sx,"sy"=sb$sy,"sum"=sb$sum,"slack"=sb$slack)
}

\keyword{models}

