% $Id: Benchmarking-package.Rd 72 2010-09-11 17:06:14Z Lars $
\name{Benchmarking-package}
\alias{Benchmarking-package}
\alias{Benchmarking}
\docType{package}

\title{ Benchmarking: Benchmark and efficiency: deterministic (DEA) and 
  stochastic frontiers (SFA) --  Estimation and plot of technologies
}

\description{

Calculates efficiency for DEA models under various technology
assumptions, allowing for slacks if wanted

Calculates super efficiency in DEA models

Bootstrap of DEA models, test in DEA models

Malmqvist index for DEA models

Estimate stochastic frontier production function (SFA)

Plot of DEA technologies for two goods; either a transformation curve
(2 outputs), isoquant (2 inputs), or a production function  (1 input og 1
output).

}

\details{
\tabular{ll}{
Package: \tab Benchmarking\cr
Type: \tab Package\cr
Version: \tab 0.10 ($Revision: 72 $)\cr
Date: \tab $Date: 2010-09-11 19:06:14 +0200 (11 sep 2010) $\cr
License: \tab Copyright\cr
}

\tabular{ll}{
\code{dea} \tab Estimate DEA input or output efficience measures,
peers, lambdas and slacks \cr

\code{slack} \tab {calculate slacks in DEA models} \cr

\code{dea.dual} \tab calculates dual weights (prices), including restrictions on weights \cr

\code{dea.merge}\tab calculate and decompose efficiency from a merger
of firms\cr

\code{sdea} \tab calculates super efficiency.\cr
\code{sfa} \tab estimate stochastic frontier production function \cr

\code{cost.opt} \tab Calculate the optimal input for given output and
prices; accept negative input and output.\cr

\code{revenue.opt} \tab Calculate the optimal output for given input
and prices; accept negative input and output.\cr

\code{profit.opt} \tab Calculate the optimal input and output for
given input and output prices.\cr

\code{dea.plot} \tab making DEA graphs under alternative
 technology assumptions.\cr

\code{dea.plot.frontier} \tab is specialized for 1 input and 1 output.\cr
\code{dea.plot.isoquant} \tab is specialized for 2 inputs.\cr
\code{dea.plot.transform} \tab is specialized for 2 outputs.\cr

\code{critValue}\tab Critical values calculated from bootstraped DEA
models. \cr

\code{typeIerror}\tab Probability of a type I error for a test in
bootstraped DEA models. \cr 
}
}

\author{
Lars Otto \email{lo@foi.dk}
}

\references{
Andersen & Petersen (1993)\cr
FEAR (Wilson 2008)\cr
}

\note{ The intterface for the methods looks very much like the
interface to the methods in the package FEAR (Wilson 2008). The
changes are that the data here are transposed compared to FEAR for in
this way to look like the way R handles data, the argument for the
options \code{RTS} and \code{ORIENTATION} can also be given as
memotechnical strings, and there are more options to control
output. The return is now not just the efficiency, but also slacks,
dual values (shadow prices), peers, and lambdas.  }

\seealso{
FEAR (Wilson 2008)
}

\examples{
# Plot of different technologies
x <- matrix(c(100,200,300,500),ncol=1,dimnames=list(LETTERS[1:4],"x"))
y <- matrix(c(75,100,300,400),ncol=1,dimnames=list(LETTERS[1:4],"y"))
dea.plot(x,y,RTS="vrs",ORIENTATION="in-out",txt=rownames(x))
dea.plot(x,y,RTS="drs",ORIENTATION="in-out",add=TRUE,lty="dashed",lwd=2)
dea.plot(x,y,RTS="crs",ORIENTATION="in-out",add=TRUE,lty="dotted")
                      
                      
dea.plot(x,y,RTS="fdh",ORIENTATION="in-out",txt=rownames(x),main="fdh")
dea.plot(x,y,RTS="irs",ORIENTATION="in-out",txt=rownames(x),main="irs")
dea.plot(x,y,RTS="irs2",ORIENTATION="in-out",txt=rownames(x),main="irs2")
# dea.plot(x,y,RTS="add",ORIENTATION="in-out",txt=rownames(x))


#  A quick frontier with 1 input and 1 output
dea.plot(x,y, main="Basic plot of frontier")


# Calculating of efficiency
dea(x,y, RTS="vrs", ORIENTATION="in")
e <- dea(x,y, RTS="vrs", ORIENTATION="in")
e
eff(e)
peers(e)
matrix(rownames(x)[peers(e)],ncol=dim(x)[1])
print(
	matrix(rownames(x)[peers(e)],ncol=dim(x)[1], 
	dimnames=list(c( paste("r",seq(1,dim(peers(e))[2]),sep="") ), rownames(x)) )
, quote=FALSE)
lambda(e)


# Calculating of super efficiency
sdea(x,y, RTS="vrs", ORIENTATION="in")
# Technology for super efficiency for firm number 3/C 
dea.plot(x[-3],y[-3],RTS="vrs",ORIENTATION="in-out",txt=LETTERS[c(1,2,4)])
points(x[3],y[3],cex=2)
text(x[3],y[3],LETTERS[3],adj=c(-.75,.75))
e3 <- dea(x,y, XREF=x[-3,], YREF=y[-3,])
eff(e3)
peers(e3)
print(matrix(rownames(x)[peers(e3)],ncol=dim(x)[1]),quote=FALSE)
lambda(e3)
e3$lambda



# Taking care of slacks
x <- matrix(c(100,200,300,500,100,600),ncol=1,
        dimnames=list(LETTERS[1:6],"x"))
y <- matrix(c(75,100,300,400,50,400),ncol=1,
        dimnames=list(LETTERS[1:6],"y"))

# Phase one
e <- dea(x,y)
print(e)
peers(e)
lambda(e)
# Phase two
sl <- slack(x,y,e)
data.frame(sl$sx,sl$sy)
peers(sl)
lambda(sl)
sl$lambda

# Slacks, two phases in one
e2 <- dea(x,y,SLACK=TRUE)
print(e2)
data.frame(eff(e2),e2$slack,e2$sx,e2$sy,lambda(e2))
peers(e2)
lambda(e2)
e2$lambda

}


\keyword{package}
\keyword{models}
\keyword{efficiency}

