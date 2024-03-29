% $Id: malmquist.Rd 247 2022-09-15 22:38:27Z X052717 $

\name{malmquist}
\alias{malmquist}

\title{Malmquist index for firms in a panel}

\description{Estimate Malmquist index for firms in a panel data set. The data set does not
    need to be balanced.}

\usage{
malmquist(X, Y, ID, TIME, RTS = "vrs", ORIENTATION = "in", SAMEREF=FALSE, 
    SLACK = FALSE, DUAL = FALSE, DIRECT = NULL, param = NULL, 
    TRANSPOSE = FALSE, FAST = TRUE, LP = FALSE, CONTROL = NULL, LPK = NULL)
}

\arguments{

  \item{X}{Inputs of firms in many periods, a (T*K) x m matrix
  of observations of K firms with m outputs (firm x input) in at the most T periods.}
  
  \item{Y}{Outputs of firms in many periods, a (T*K) x n matrix
  of observations of K0 firms with n outputs (firm x input) in at the most T periods.}
  
  \item{ID}{Identifier for the firms in rows of \code{X} and \code{Y}.}

  \item{TIME}{Array with period number for each row in the input maxtrix \code{X} and output matrix\code{Y}}
  
  \item{RTS}{Returns to scale assumption as in \code{\link{dea}}.}
    
  \item{ORIENTATION}{Input efficiency "in" (1), output
    efficiency "out" (2), and graph efficiency "graph" (3) as in \code{\link{dea}}.}

  \item{SAMEREF}{Use the same units for reference technology when comparing 
  two periods. This is not restricted to same units in several timpe periods, 
  but only to pairwise periods comparisons for Malmquist. Default is to use 
  available and possible differnt units in pairwise periods.}

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

\details{\code{Malmquist} uses \code{\link{malmq}} for the calculations of the 
   necessary efficiencies, and the returned indices are as in \code{\link{malmq}}.
   The data must be a long data set with regards to \code{TIME} and \code{ID}; se the example below.
   
   Note that the calculated index are index comparing a period and the previous period. To compare the 
   development over time the indices must be turned into a chain index as shown in the example below.
}


\value{
   \item{m}{Malmquist indicies, an array of length T*K in the order of \code{ID} and \code{TIME}, 
    i.e. the order of the rows of \code{X}.}
   \item{tc}{Technical change indices, an array of length T*K.}
   \item{ec}{Efficiency indices, an array of length T*K.}
   \item{id}{Index for firms as \code{ID}}
   \item{time}{Index for time as \code{TIME}}
   \item{e00}{The efficiencies for period 0 with reference technology from period 0.}
   \item{e10}{The efficiencies for period 1 with reference technology from period 0.}
   \item{e11}{The efficiencies for period 1 with reference technology from period 1.}
   \item{e01}{The efficiencies for period 0 with reference technology from period 1.}
}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
DEA, SFA, and R}; Springer 2011 }

\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}

\note{
   The lagged values \code{e11} are not necessary equal to values of \code{e00}
   as the reference technology for the two periods could be generated by
   different units, if the units in different time periods are not the same.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dea}}, \code{\link{malmq}}}

\examples{
 x0 <- matrix(c(10, 28, 30, 60),ncol=1)
 y0 <- matrix(c(5, 7, 10, 15),ncol=1)
 x1 <- matrix(c(12, 26, 16, 60 ),ncol=1)
 y1 <- matrix(c(6, 8, 9, 15 ),ncol=1)
 x2 <- matrix(c(13, 26, 15, 60 ),ncol=1)
 y2 <- matrix(c(7, 9, 10, 15 ),ncol=1)
 
 dea.plot(x0, y0, RTS="vrs", txt=TRUE)
 dea.plot(x1, y1, RTS="vrs", add=TRUE, col="red")
 dea.plot(x2, y2, RTS="vrs", add=TRUE, col="blue")
 points(x1, y1, col="red", pch=16)
 # points(x2, y2, col="blue", pch=17)
 text(x1, y1, 1:dim(x1)[1], col="red", adj=-1)
 text(x2, y2, 1:dim(x1)[1], col="blue", adj=-1)
 legend("bottomright", legend=c("Period 0", "Period 1", "Period 2"),
    col=c("black", "red", "blue"), lty=1, pch=c(1,16, 17),  bty="n")
 
 X <- rbind(x0, x1, x2)
 Y <- rbind(y0, y1, y2)
 # Make ID and TIME variables one way or another
 ID <- rep(1:dim(x1)[1], 3)
 # TIME <- c(rep(0,dim(x1)[1]), rep(1,dim(x1)[1]), rep(2,dim(x1)[1]))
 TIME <- gl(3, dim(x1)[1], labels=0:2)
 # This is how the data for Malmquist must look like
 data.frame(TIME, ID, X, Y)
 mq <- malmquist(X,Y, ID, TIME=TIME) 
 data.frame(TIME, ID, X, Y, mq$e00, mq$e01, mq$e10, mq$e11, mq$m, mq$tc)[order(ID, TIME),]
 
 # How to make the Malmquist indices to a chain index
 # Make data.frame with indices
 DM <- data.frame(TIME, ID, m=mq$m, tc=mq$tc, ec=mq$ec)
 # Set missing index for first period to 1, the base
 DM[DM$TIME==0, c("m","tc", "ec")] <- 1
 # Make chain index of the individual indices
 AD <- aggregate(cbind(m=DM$m), by=list(ID=DM$ID), cumprod)
 # Compare chain index to original index
 data.frame(ID, TIME, m=c(AD$m), DM$m)
 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ efficiencies }
\keyword{ Malmquist }
\keyword{ panel }
