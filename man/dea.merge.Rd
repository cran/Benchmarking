% $Id: dea.merge.Rd 72 2010-09-11 17:06:14Z Lars $
\name{dea.merge}
\alias{dea.merge}

\title{
Efficiency results after mergers and their decompositions
}

\description{ Calculate and decompose efficiency from mergers of
similar firms, horisontal integration.  }

\usage{

dea.merge(X, Y, M, RTS = "vrs", ORIENTATION = "in",
             XREF = NULL, YREF = NULL, FRONT.IDX = NULL, TRANSPOSE = FALSE,
             LP = FALSE, CONTROL = NULL, LPK = NULL, ...)
}

\arguments{ Most of the arguments corresponds to the argument in
\code{dea}, K firms, m inputs, and n outputs.

  \item{X}{ K times m matrix as in \code{dea} }

  \item{Y}{ K times n matrix as in \code{dea} }

  \item{M}{Kg times K matrix where each row defines a merger by the
  firms ((colloms) included; matrix as returned from method
  \code{make.marge}}

  \item{RTS}{as in \code{dea} }

  \item{ORIENTATION}{ as in \code{dea} }

  \item{XREF}{ as in \code{dea} }

  \item{YREF}{ as in \code{dea} }

  \item{FRONT.IDX}{ as in \code{dea} }

  \item{TRANSPOSE}{ as in \code{dea} }

  \item{LP}{ as in \code{dea} }

  \item{CONTROL}{ as in \code{dea} }

  \item{LPK}{ as in \code{dea} }

  \item{\dots}{ as in \code{dea} }

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{

\item{Eff}{overall efficiencies of mergers, Kg vector}

\item{Estar}{adjusted overall efficiencies of mergers after the
removal of individual learning, Kg vector}

\item{learning}{= individual learning effects me (row) vectors, K*
(row) vectors}

\item{harmony}{harmony effects, Kg vector}

\item{size}{size effects, K* vector0}

}

\references{
PB & LO, chapter XX
}

\author{
Peter Bogetoft and Lars Otto \email{lo@foi.dk}
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
