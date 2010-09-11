% $Id: dea.plot.Rd 72 2010-09-11 17:06:14Z Lars $
\name{dea.plot}
\alias{dea.plot}
\alias{dea.plot.frontier}
\alias{dea.plot.isoquant}
\alias{dea.plot.transform}


\title{Plot of DEA technologies}

\description{Draw a graph of DEA technologies given various
assumptions on scale of returns. Only for two goods, either isoquant
(2 inputs), transformationscurve (2 outputs), or a production function
(1 input and 1 output).  If the number of good is larger than 2 then
aggragation occur, either simple or weighted.}

\usage{

dea.plot(x, y, RTS="vrs", ORIENTATION="in-out", txt=NULL, add=FALSE, 
            wx=NULL, wy=NULL, TRANSPOSE = FALSE, fex=1,
            ..., RANGE=FALSE, xlim, ylim, xlab, ylab)

dea.plot.frontier(x, y, RTS="vrs",...)

dea.plot.isoquant(x1, x2, RTS="vrs",...)

dea.plot.transform(y1, y2, RTS="vrs",...)

}


\arguments{

  \item{x}{The good for the first axis}

  \item{y}{The good for the second axis}
  \item{x1, y1}{The good for the first axis}
  \item{x2, y2}{The good for the second axis}

  \item{RTS}{Assumptions for returns to scale: "fdh" (0), "vrs" (1),
      "drs" (2), "crs" (3), "irs" (4), "irs2" (5) (irs without
      convexity).  To be added later "add" (6).  Numbers in
      parenthesis can also be used as values for \code{RTS}}

  \item{ORIENTATION}{Input-output graph 1 input and 1 output "in-out"
      (0), graph of 2 inputs "in" (1), and graph of 2 outputs "out"
      (2). }



  \item{txt}{\code{txt} is an array to label the observations}

  \item{add}{For \code{add=T} the technology is drawn on top of an 
    existing graph, and for the defautlt \code{add=F} a new graph is made.}

  \item{wx}{Weight to aggregate for the first axis if there are more
  than 1 good in the }
  
  \item{wy}{Weights to aggregater for the second axis if there is more
  than 1 good }

  \item{TRANSPOSE}{Only relevant for more than 1 input for each axis,
  see \code{dea} for a description}

  \item{\dots}{Usual options for the methods \code{plot}, \code{lines}, 
      and \code{abline} etc.}

  \item{fex}{Relative size for text on observations, corresponds to 
      \code{cex}, but only changes this size.}

   \item{RANGE}{A logical variable, if \code{RANGE=TRUE} the the
   limits for the graph is the range of the variables; zero is always
   included. Default is \code{RANGE=FALSE} when the range is from zero
   to the max values.  Relevant if some values are negative.}

  \item{xlim}{Possible limits \code{c(x1,x2)} for the first axis}

  \item{ylim}{Possible limits  \code{c(y1,y2)} for the second axis}

  \item{xlab}{Possible label for the x-axis}

  \item{ylab}{Possible label for the y-axis}

}


\details{ The method \code{dea.plot} is the general plottting method,
  the 3 other are specialized version for frontiers, 1 input and 1
  output, isoquant curves, 2 inputs for given outputs, and
  transformation curves, 2 outputs for given inputs.  }


\value{No return, using the original graphing system.
}

\references{
Uses classical R graphics.  Should be implementet in grid.
}

\author{
Lars Otto 
}

\note{
No notes
}


\seealso{
\code{plot}
}

\examples{

x <- matrix(c(100,200,300,500,600,100),ncol=1)
y <- matrix(c(75,100,300,400,400,50),ncol=1)

dea.plot(x,y,RTS="vrs",ORIENTATION=3,txt=LETTERS[1:length(x)])
dea.plot(x,y,RTS="crs",ORIENTATION=3,add=TRUE,lty="dashed")

dea.plot.frontier(x,y,txt=1:dim(x)[1])

n <- 10
x <- matrix(1:n,,1)
y <- matrix(x^(1.6) + abs(rnorm(n)),,1)
dea.plot.frontier(x,y,RTS="irs",txt=1:n)
dea.plot.frontier(x,y,RTS="irs2",add=TRUE,lty="dotted")

id <- sample(1:n,30,replace=TRUE)
dea.plot(x[id],y[id],RTS="irs",ORIENTATION=3)
dea.plot.frontier(x[id],y[id],RTS="irs")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{plot}
