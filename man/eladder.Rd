% $Id: eladder.Rd 86 2010-11-11 09:00:13Z Lars $
\name{eladder}
\alias{eladder}
\alias{eladder.plot}

\title{ Efficiency ladder for a single firm }

\description{ How the efficiency changes as the most influential peer
is removed sequentially one at a time }

\usage{
eladder(n, X, Y, RTS = "vrs", ORIENTATION = "in")
eladder.plot(elad, peer, TRIM = NULL)
}


\arguments{

  \item{n}{The number of the firm where the ladder is calculated }

  \item{X}{Inputs of firms to be evaluated, a K x m matrix of
      observations of K firms with m inputs (firm x input).  In case
      \code{TRANSPOSE=TRUE} the input matrix is transposed to input x
      firm. }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix of
      observations of K firms with n outputs (firm x input).  In case
      \code{TRANSPOSE=TRUE} the output matrix is transposed to output
      x firm. }

  \item{RTS}{ Text string or a number defining the underlying DEA
      technology / returns to scale assumption, se the possible values
      for \link{dea}. }

  \item{ORIENTATION}{Input efficiency "in" (1), output efficiency
      "out" (2), and graph efficiency "graph" (3). For use with
      \code{DIRECT}, an additional option is "in-out" (0).}

  \item{elad}{The sequence of efficiencies returned from \code{eladder}. }

  \item{peer}{The sequence of peers returned from \code{eladder}. }

  \item{TRIM}{The number of characters for the name of the peers on
      the axis in the plot. }

}


\details{ The function \code{eladder} calculates how the efficiency
   for a firm changes when the most influential peer is removed
   sequentially one at a time. Somewhere in the sequence the firm
   becomes efficient and are itself removed from the set of firms
   generating the technology (or the only firm left) and thereafter
   the efficiencies are super-efficiencies.}


\value{The object returned from \code{eladder} is a list with components
\item{eff}{The sequence of efficiencies}
\item{peer}{The peers as index in the \code{X} rows.}
} 


\note{When the number of firms is large then the number of influential
   peers will also be large and the names or numbers of the peers on
   the x-axis might be squeeze together and be illegible. In this case
   restrict the number of influential peers to be removed.

   The efficiency step ladder is discussed in Essay 4 of Dag Fjeld
   Edvardsens's phd thesis from 2004.}


\references{ Dag Fjeld Edvardsen; \emph{Four Essays on the Measurement
   of Productive Efficiency}; University of Gothenburg 2004;
   \url{http://hdl.handle.net/2077/2923}}


\author{ Peter Bogetoft and Lars Otto \email{lo@foi.dk} }


\examples{
data(charnes1981)
attach(charnes1981)
x <- cbind(x1,x2,x3,x4,x5)
y <- cbind(y1,y2,y3)

# Choose the firm for analysis, we choose 'Tacoma'
n <- which(name=="Tacoma")[1]

el <- eladder(n, x, y, RTS="crs")
eladder.plot(el$eff, el$peer)
eladder.plot(el$eff, name[el$peer])

# Restrict to 20 most influential peers for 'Tacoma'
eladder.plot(el$eff[1:20], name[el$peer][1:20])

# Truncate the names of the peers and put a title on top
eladder.plot(el$eff, name[el$peer], TRIM=5)
title("Eladder")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
