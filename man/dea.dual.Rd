% $Id: dea.dual.Rd 72 2010-09-11 17:06:14Z Lars $
\name{dea.dual}
\alias{dea.dual}
% \alias{dea.dual.mat}
\encoding{UTF-8}

\title{
Calculates dual DEA models
}

\description{
Calculates dual DEA models including restrictions on the ratios, acceptance regions
}

\usage{
dea.dual(X, Y, RTS = "vrs", ORIENTATION = "in", XREF = NULL, YREF = NULL, 
           FRONT.IDX = NULL, DUAL = NULL, 
           TRANSPOSE = FALSE, LP = FALSE, CONTROL=NULL, LPK=NULL, ...)


% dea.dual.mat(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL,YREF=NULL,
%         FRONT.IDX=NULL, DUAL=NULL, # SLACK=FALSE, delta=0, DELTA=FALSE, 
%         TRANSPOSE=FALSE, LP=FALSE, ...)

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

  \item{ORIENTATION}{Input efficiency "in" (1), output efficiency
      "out" (2), and graph efficiency "graph" (3) (not yet
      implemented). }


  \item{XREF}{input determining the technology, defaults to \code{X}}
  \item{YREF}{output determining the technology, defaults to \code{Y}}
  \item{FRONT.IDX}{Index for units determining the technology}


  \item{DUAL}{

 Matrix of order "number of input plus number of outputs minus 2" times
 2.  First column is the lower bound and the second column is the
 upper bound for the restrictions on the multiplier ratios inputs and
 outputs. The rations are relative to the first input and first output
 implying that there is no restriction for neither the first input nor
 the first output so the number of restrictions is two less the total
 number of inputs and outputs.


    DUALIN og DUALOUT eller blot DUAL.  Her er valgt DUAL saa DUAL
    skal vaere en matrix saaledes at der en raekke i DUAL for hvert
    input og hvert output, bortset fra det foerste input og det
    foerste output, dvs m-1+n-1=m+n-2 raekker.  Der skal vaere 2
    soejler, foerste soejle er den nedre graense, og anden soejle er
    den oevre graense.

}

   \item{TRANSPOSE}{Input and output matrices are treated as goods
      times firms as defaults.  If \code{TRANSPOSE=FALSE} then marices
      are firms times goods as is standard in R for statistical
      models}

   \item{LP}{Only for debugging. If \code{LP=TRUE} then input and
      output for the LP program are written to standard output for
      each unit.}

  \item{CONTROL}{Possible controls to lpSolveAPI, se the documentatin
  for that package.}

   \item{LPK}{When \code{LPK=k} then a mps file is written for firm
   \code{k}; it can be used as input to an alternative LP solver just
   to check the our results.}

   \item{\dots}{Optional parameters for something.}

}

\details{ Solved as an LP program using the package lpSolveAPI. The
method \code{dea.dual.dea} calls the method \code{dea} with the option
\code{DUAL=TRUE}.  }

\value{

\item{eff}{ The efficiencies }

\item{objval}{ The objcetive value as returned from, normally just as
eff }

\item{RTS}{ The return to scale assumption, just the option \code{RTS}
in the call }

\item{ORIENTATION}{ The efficiency orientation as in the call }

\item{TRANSPOSE}{ As in the call }

\item{u}{Dual values, prices, for inputs }

\item{v}{Dual values, prices, for outputs }

\item{gamma}{ The values of gamma, the shadow price(s) for returns to
scale }

\item{sol}{Solution of all variables as one component, sol=c(u,v,gamma).}

}


\references{
See \code{dea}
}

\author{
Lars Otto \email{lo@foi.dk}
}

\note{ 

Note that the dual values are not unique for extreme points in the
technology set, i.e. for firms where efficiency is 1 and there is no
slack. In this case the value of the calculated dual variable can
depend on the order of the complete efficient firms.

}



\seealso{
\code{dea}
}

\examples{

x <- matrix(c(2,5 , 1,2 , 2,2 , 3,2 , 3,1 , 4,1), ncol=2,byrow=TRUE)
y <- matrix(1,nrow=dim(x)[1])
dea.plot.isoquant(x[,1],x[,2],txt=1:dim(x)[1])
for (i in 1:dim(x)[1]) lines(c(0,x[i,1]), c(0,x[i,2]),lty="dotted")

e <- dea(x,y,RTS="crs",SLACK=TRUE)
ed <- dea.dual(x,y,RTS="crs")
# print(cbind("e"=e$eff,"ed"=ed$eff, peers(e), lambda(e), e$sx, e$sy, ed$u, ed$v),
#     digits=3)
# 
# (dual <- matrix(c(.5, 2.5), nrow=dim(x)[2]-1, ncol=2, byrow=TRUE))
# er <- dea.dual(x,y,RTS="crs", DUAL=dual)
# print(cbind("e"=e$eff,"ar"=er$eff, lambda(e), e$sx, e$sy, er$u, 
#     "ratio"=er$u[,2]/er$u[,1],er$v),digits=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
