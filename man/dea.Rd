% $Id: dea.Rd 72 2010-09-11 17:06:14Z Lars $
\name{dea}
\alias{dea}
% \alias{dea.mat}
\alias{efficiency}
\alias{effi}
\alias{eff}
\alias{sdea}
\alias{slack}
\alias{print.slack}
\alias{summary.slack}
\alias{print.Farrell}
\alias{summary.Farrell}
\alias{peers}
\alias{get.number.peers}
\alias{peerslambda}
\alias{print.peers}
\alias{lambda}
\alias{lambda.print}
\encoding{UTF-8}

\title{DEA efficiency}

\description{Estimate a DEA frontier and calculate input or output
efficience measures, peers, lambdas and slacks}

\usage{
dea(X, Y, RTS = "vrs", ORIENTATION = "in", XREF = NULL, YREF = NULL,
    FRONT.IDX = NULL, SLACK = FALSE, DUAL=TRUE,
    TRANSPOSE = FALSE, FAST = FALSE, LP = FALSE, CONTROL=NULL, 
    LPK=NULL, ...)

efficiency(object)
eff(object)

print.Farrell(x, digits=4, ...) 

summary.Farrell(object, digits=4, ...)

slack(X, Y, e, XREF=NULL, YREF=NULL, FRONT.IDX=NULL, LP=FALSE)

peers(object,  NAMES=FALSE)
get.number.peers(object)
peerslambda(object)
print.peers(x, ...)
lambda(object, KEEPREF=FALSE, ...)
lambda.print(x, KEEPREF=FALSE, ...)

sdea(X, Y, RTS="vrs", ORIENTATION="in", TRANSPOSE=FALSE, LP=FALSE, ...)  

}



%- maybe also 'usage' for other objects documented here.
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
      "out" (2), and graph efficiency "graph" (3). }

  \item{XREF}{input determining the technology, defaults to \code{X}}
  \item{YREF}{output determining the technology, defaults to \code{Y}}
  \item{FRONT.IDX}{Index for units determining the technology}

  \item{SLACK}{Calculate slack in a phase II calculation; the
  precision for clalculating slacks for orientation graph is low.}

  \item{DUAL}{Calculate dual variables, i.e. shadow prices; not
  calculated for orientation graph as that is not an LP problem.}

   \item{TRANSPOSE}{Input and output matrices are treated as goods
      times firms as defaults.  If \code{TRANSPOSE=FALSE} then marices
      are firms times goods as is standard in R for statistical
      models}

   \item{LP}{Only for debugging. If \code{LP=TRUE} then input and
      output for the LP program are written to standard output for
      each unit.}

   \item{FAST}{Only calculate efficiencies and just return them as a
   vector, i.e. no lambda or other outout. The return when using FAST
   cannot be used as input for \code{slack}.}

  \item{CONTROL}{Possible controls to lpSolveAPI, se the documentatin
  for that package.}

   \item{\dots}{Optional parameters for something.}

   \item{object, x}{A Farrel object, the object returned from \code{dea}}

   \item{KEEPREF}{If \code{TRUE} then all firms are kept as reference
   firms even though they have all zero weights (lambda); might come
   handy if one needs to calculate \code{X x lambda} such that the
   firms in \code{X} and \code{lambda} agree. If \code{FALSE}, the
   default, then only weight for the peers are in the matrix
   \code{lambda}. }

   \item{NAMES}{Use names if any in the return for \code{peers}.}

   \item{e}{Object of type Farrell returned from \code{dea}.}

   \item{digits}{Digits in printes output, handled by format in print.}

   \item{LPK}{When \code{LPK=k} then a mps file is written for firm
   \code{k}; it can be used as input to an alternative LP solver just
   to check the our results.}

}

\details{ The efficiency in \code{dea} is calculated by the LP method
in the package lpSolveAPI.  Slacks can be calculated either in the
call of \code{dea} using the option \code{SLACK=TRUE} or in a
following call to the method \code{slack}

The methods \code{sdea} calculate super-effciency and returns the same
type of object as \code{dea}.

The \code{peers} returns a matrix with numbers for the peers for each
firm, for firms with efficiency 1 the peers are the firm itself.  If
there is slack for a firm with efficiency 1 then it can be questioned
if it is reasonbable to consicer the firm as a real peer, and this can
be corrected with a call to slack, either directly or by option SLACK.

The \code{get.number.peers} returns the number of firms that a peer is
a peer of.

% The method \code{dea.mat} uses \code{lpSolve} as the LP engine instead
% of the \code{lpSolveAOI} LP engine. The result is that \code{dea.mat}
% is much slower than \code{dea} except for data with very small number
% of firms, like less than 10 firms.

The graph orientated efficiency is calculated by bisection between
feasible and not feasible values of G. The precision in the result is
less than for the other orientations.

}


\value{ The results are returned in a Farrell object with the
following components.  The method \code{slack} only returns the last three
componnets in the list relevant for slacks.

\item{eff}{The efficiencies }

\item{lambda}{The lambdas, i.e. the weight of the peers for each firm}

\item{objval}{The objcetive value as returned from, normally just as
eff, for for \code{slack} the sum of the slacks. }

\item{RTS}{The return to scale assumption, just the option \code{RTS}
in the call}

\item{ORIENTATION}{The efficiency orientation as in the call}

\item{TRANSPOSE}{As in the call}

\item{slack}{A vector with sums of the slacks for each firm.}

\item{sx}{A matrix for input slacks for each firm, only calculated if
the option \code{SLACK} is \code{TRUE} or returned from the method
\code{slack}}

\item{sy}{A matrix for output slack, see \code{sx}}

\item{vx}{Dual variable for input, only calculated if \code{DUAL} is
\code{TRUE}.}

\item{uy}{Dual variable for output, only calculated if \code{DUAL} is
\code{TRUE}.}

}


\references{
Bogetoft and Otto 
}

\author{
Lars Otto \email{lo@foi.dk}
}

\note{ 

Calculate slacks for super efficiency should be done by using the
option \code{SLACK=TRUE} in the call of the method \code{sdea}.  If
the two phases are done in two steps as first a call to \code{sdea}
and then a call to \code{slacks} the user must make sure to set the
reference technology to the one coresponding for super-efficiency in
the call to \code{slack} and this would need a loop with calls to
\code{slack}.

When there is slack, and slack is not taken into consideration, then
the peers for a firm with slack might depend on the order of firms in
the data set; this is a property of the LP algorithm used to solve the
problem.

}




\seealso{
FEAR 
}

\examples{

x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)
dea.plot.frontier(x,y,txt=1:dim(x)[1])

e <- dea(x,y)
effi(e)
print(e)
summary(e)
lambda(e)
peers(e)

# Input savings potential for each firm
(1-eff(e)) * x
(1-e$eff) * x

# calculate slacks
(sl <- slack(x,y,e))
data.frame(e$eff,sl$slack,sl$sx,sl$sy)

el <- dea(x,y,SLACK=TRUE)
data.frame(e$eff,el$eff,el$slack,el$sx,el$sy)

# Fully efficient units, eff==1 and no slack
which(eff(e) == 1 & !sl$slack)


# Super efficiency
se <- sdea(x,y)
data.frame(eff(e),eff(se), lambda(se))
rngx <- range(x) * c(0,1.2)
rngy <- range(y) * c(0,1.2)
n <- 3
dea.plot.frontier(x[-n],y[-n],txt=(1:4)[-n],lty=1, xlim=rngx, ylim=rngy)
points(x[n],y[n],cex=1.25,pch=16)
text(x[n],y[n],n,adj=c(-.75,.75))




}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ efficiency }% __ONLY ONE__ keyword per line
