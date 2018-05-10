% $Id: malmq.Rd 181 2018-05-06 14:58:02Z lao $

\name{malmq}
\alias{malmq}

\title{Malmquist index}

\description{
	Estimates Malmquist indices for productivity and its decomposition beteween two
	periods. The units in the two periods does not have to be exactly the same, but
	the Malmquist index is only calculated for units present in both periods.
}

\usage{
malmq(X0, Y0, ID0 = NULL, X1, Y1, ID1 = NULL, RTS = "vrs", ORIENTATION = "in", 
	SLACK = FALSE, DUAL = FALSE, DIRECT = NULL, param = NULL, TRANSPOSE = FALSE,
	FAST = TRUE, LP = FALSE, CONTROL = NULL, LPK = NULL)
}

\arguments{
  \item{X0}{Inputs of firms in period 0, a K0 x m matrix
    of observations of K0 firms with m inputs (firm x input). }
    
  \item{Y0}{Outputs of firms in period 0, a K0 x n matrix
  of observations of K0 firms with n outputs (firm x input). }
  
  \item{ID0}{Index for firms in period 0; could be numbers or labels. Length K0.}
  
  \item{X1}{Inputs of firms in period 1, a K1 x m matrix
    of observations of K1 firms with m inputs (firm x input).}
    
  \item{Y1}{Outputs of firms in period 1, a K1 x n matrix
  of observations of K1 firms with n outputs (firm x input). }
  
  \item{ID1}{Index for firms in period 0; could be numbers or labels. Length K0.}
  
  \item{RTS}{Returns to scale assumption as in \code{\link{dea}}.}
  
  \item{ORIENTATION}{Input efficiency "in" (1), output
      efficiency "out" (2), and graph efficiency "graph" (3) as in \code{\link{dea}}.}
      
  \item{SLACK}{See \code{\link{dea}}.}
  \item{DUAL}{See \code{\link{dea}}.}
\item{DIRECT}{See \code{\link{dea}}.}
\item{param}{See \code{\link{dea}}.}
\item{TRANSPOSE}{See \code{\link{dea}}.}
\item{FAST}{See \code{\link{dea}}.}
\item{LP}{See \code{\link{dea}}.}
\item{CONTROL}{See \code{\link{dea}}.}
\item{LPK}{See \code{\link{dea}}.}
  
}

\details{
See \dots\ for details
}

\value{
	\item{m}{Malmquist index for productivity.}
	\item{tc}{Index for technoligy change.}
	\item{ec}{Index for efficiency change.}
	\item{mq}{Malmquist index for productivity; same as \code{m}.}
	\item{id}{Index for firms present in both period 0 and period 1.}
	\item{id0}{Index for firms in period 0 that are also in period 1.}
	\item{id1}{Index for firms in period 1 that are also in period 0.}
	\item{e00}{The efficiencies for period 0 with reference technoligy from period 0.}
	\item{e10}{The efficiencies for period 1 with reference technoligy from period 0.}
	\item{e11}{The efficiencies for period 1 with reference technoligy from period 1.}
	\item{e01}{The efficiencies for period 0 with reference technoligy from period 1.}
}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
DEA, SFA, and R}; Springer 2011 }

\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{
   The calculations of efficiencies are only done for units present in both periods.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dea}}}

\examples{
   x0 <- matrix(c(10, 28, 30, 60),ncol=1)
   y0 <- matrix(c(5, 7, 10, 15),ncol=1)
   x1 <- matrix(c(12, 26, 16, 60 ),ncol=1)
   y1 <- matrix(c(6, 8, 9, 15 ),ncol=1)

   dea.plot(x0, y0, RTS="vrs", txt=TRUE)
   dea.plot(x1, y1, RTS="vrs", add=TRUE, col="red")
   points(x1, y1, col="red", pch=16)
   text(x1, y1, 1:dim(x1)[1], col="red", adj=-1)

   m <- malmq(x0,y0,,x1,y1,,RTS="vrs")
   print("Malmquist index for change in productivity, technoligy change:")
   print(m$mq)
   print("Index for change of frontier:")
   print(m$tc)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ efficiencies }
\keyword{ Malmquist }
