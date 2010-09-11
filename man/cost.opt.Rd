% $Id: cost.opt.Rd 72 2010-09-11 17:06:14Z Lars $
\name{cost.opt}
\alias{cost.opt}
\alias{print.cost.opt}
\alias{summary.cost.opt}

\alias{revenue.opt}
\alias{print.revenue.opt}
\alias{summary.revenue.opt}

\alias{profit.opt}
\alias{print.profit.opt}
\alias{summary.profit.opt}

\title{DEA Optimal Cost}

\description{Estimate the optimal input and/or output for
optimal/minimum/maximum cost, revenue or profit for a DEA technology
set}

\usage{
cost.opt(XREF, YREF, W, YOBS=NULL, RTS="vrs", TRANSPOSE=FALSE, 
                     LP=FALSE, LPK = NULL)  

revenue.opt(XREF, YREF, P, XOBS = NULL, RTS = "vrs", TRANSPOSE = FALSE, 
              LP = FALSE, LPK = NULL)

profit.opt(XREF, YREF, W, P, RTS = "vrs", TRANSPOSE = FALSE, 
                 LP = FALSE, LPK = NULL)

}

\arguments{

  \item{}{Input and output matrices are in the same form as for the
  method \code{dea}.}

  \item{XREF}{input determining technology where efficiency is to be
  calculated, a K x m matrix of observations of K firms with m inputs
  (firm x input).  In case \code{TRABSPOSE=TRUE} the input matrix is
  transposed as input x firm. }

  \item{YREF}{output determinig technology where efficiency is to be
  calculated, a K x n matrix of observations of K firms with n outputs
  (firm x input).  In case \code{TRABSPOSE=TRUE} the output matrix is
  transposed as output x firm. }

  \item{W}{input prices as a matrix. Either same prices for all firms
  or individual prices for all firms; i.e. etiher a 1 x m or K x m
  matrix for K firms and m inputs }

 \item{P}{output prices as a matrix. Either same prices for all firms
  or individual prices for all firms; i.e. etiher a 1 x n or K x n
  matrix for K firms and n outputs }

  \item{XOBS}{the input to where optimal output is to be
  calculated. Deaults to \code{XREF}.  Same form as \code{XREF}}

  \item{YOBS}{the output to where optimal input is to be
  calculated. Deaults to \code{YREF}.  Same form as \code{YREF}}

  \item{RTS}{

  A text string or a number to state the returns to scale
  assumption.  Only numbers implementet as of this writing.

\tabular{rcl}{
0 \tab "fdh" \tab Free disposability hull, no convexity assumption \cr
1 \tab "vrs" \tab Variable returns to scale, convex and free disposability \cr
2 \tab "drs" \tab Decreasing returns to scale, convex, down-scaling and free disposability \cr
3 \tab "crs" \tab constant returns to scale and free disposability \cr
4 \tab "irs" \tab increasing returns to scale, up-scaling, but not down-scaling \cr
5 \tab "add" \tab additivity,  scaling up and down, but only with integers
}

}

  \item{TRANSPOSE}{Input and output matrices are treated as goods
        times units as defaults.  If \code{TRANSPOSE=FALSE} then
        marices are units times goods as in standard in R for
        statistical models}

   \item{LP}{Only for debugging. If \code{LP=TRUE} then input and
      output for the LP program are written to standard output for
      each unit.}

   \item{LPK}{When \code{LPK=k} then a mps file is written for firm
   \code{k}; it can be used as input to an alternative LP solver just
   to check the our results.}

}

\details{ The optimization is done by the R method in the package
  \code{lpSolve}.  The methods \code{print} and \code{summary} is
  working for \code{cost.opt}

} 

\value{ The values returned is the optimal input. When saved in an
object the following is available:

\item{xopt}{The optimal input, returned as a matrix.}

\item{cost}{The optimal cost.}

\item{lambda}{The weights that determines the technology.}

}

\references{
PB & LO
}

\author{Lars Otto \email{lo@foi.dk}}

\note{

The index for peer units can be returned by the method \code{peers}
and the weights are returned in \code{lambda}.  Note that the peers
now are the firms for the optimal input allocation, not just the
technical efficient firms.

}


\seealso{
FEAR
}

\examples{

x <- matrix(c(2,12, 2,8, 5,5, 10,4, 10,6, 3,13), ncol=2, byrow=TRUE)
y <- matrix(1,nrow=dim(x)[1],ncol=1)
w <- matrix(c(1.5, 1),ncol=2)

txt <- LETTERS[1:dim(x)[1]]
dea.plot(x[,1],x[,2], ORIENTATION="in",  cex=1.25)
text(x[,1],x[,2],txt,adj=c(-.7,-.2),cex=1.25)

te <- dea(x,y,RTS="vrs")
xopt <- cost.opt(x,y,w,RTS=1)
cobs <- x \%*\% t(w)
copt <- xopt$x \%*\% t(w)
ce <- copt/cobs
ae <- ce/te$eff
data.frame("ce"=ce,"te"=te$eff,"ae"=ae)
cbind(ce,eff(te),ae)
print(cbind("ce"=c(ce),"te"=te$eff,"ae"=c(ae)),digits=2)

# isocost line in the technology plot
abline(a=copt[1]/w[2], b=-w[1]/w[2], lty="dashed")
abline(a=0, b=x[2,2]/x[2,1], lty="dotted")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{efficiency}
\keyword{DEA}
