% $Id: dea.boot.Rd 96 2010-11-27 22:47:25Z Lars $
\name{dea.boot}
\alias{dea.boot}


\title{ Bootstrap DEA models }

\description{ A wrapper for \code{boot.sw89} from the package
   \pkg{FEAR} such that the returns are bootstrap of Farrell
   efficiencies.  All bootstrap is done by in \code{boot.sw89}. }


\usage{
dea.boot(X, Y, NREP = 50, EFF = NULL, RTS = "vrs", ORIENTATION = "in", 
         alpha = 0.05, XREF = NULL, YREF = NULL, EREF = NULL)
}


\arguments{

  \item{X}{Inputs of firms to be evaluated, a K x m matrix of
      observations of K firms with m inputs (firm x input) }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix of
      observations of K firms with n outputs (firm x input). }

  \item{NREP}{Number of bootstrap replicats}

  \item{EFF}{Efficiencies for (X,Y) relative to the technology
      generated from (XREF,YREF).}

  \item{RTS}{The returns to scale assumptions as in \code{link{dea}},
      only works for "vrs", "drs", and "crs".}

  \item{ORIENTATION}{Input efficiency "in" (1), output
      efficiency "out" (2), and graph efficiency "graph" (3). }

  \item{alpha}{One minus the size of the confidence interval for the
      bias corrected efficiencies}

  \item{XREF}{Inputs of the firms determining the technology, defaults
      to \code{X}. }

  \item{YREF}{Outputs of the firms determining the technology,
      defaults to \code{Y}. }

  \item{EREF}{Efficiencies for the firms in XREF, YREF.}

}


\details{ The details are lightly explained in Bogetoft and Otto
   (2011) chapter 6.  

   The bootstrap at the moment does not work for any kind of
   directional efficiency.  

   The returned confidence interval is for the bias corrected
   efficiencies; to get confidence interval for the uncorrected
   efficiencies add the bias to both upper and lower values for the
   interval.}



\value{ The returned values are the values from \code{boot.sw89}
   recalculated recoded to suit Farrell efficiencies in Benchmarking.  

  \item{bias}{An array of bootstrap bias estimates for the firms in X,Y}

  \item{var}{An array of bootstrap variance estimates for the firms in
      X,Y}

  \item{conf.int}{\code{K x 2} matrix with confidence interval for the
      estimated efficiencies}

  \item{ehat}{Efficiencies}

  \item{ehat.bc}{Bias-corrected efficiencies}

%  \item{bias.flag}{ TRUE if a bias correction is negative.}

  \item{boot}{ The replicat bootstrap estimates of the efficiencies, a
      K times NREP matrix}.

}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
   DEA, SFA, and R}; Springer 2011

   Paul W. Wilson (2008), \dQuote{FEAR 1.0: A Software Package for
   Frontier Efficiency Analysis with R,} \emph{Socio-Economic Planning
   Sciences} 42, 247--254 }


\author{ Peter Bogetoft and Lars Otto \email{lo@foi.dk}}


\note{ \emph{The function \code{dea.boot} only works if the package
   \pkg{FEAR} is installed}; it does not need to be loaded.

   The \code{boot.sw89} is suited for Shephard distance functions and
   the present present function is a wrapper to recalculate the
   results suitable for Farrell efficiencies.  }


\seealso{ The documentation for \code{boot.sw98} in the package
   \pkg{FEAR}. }


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)

e <- dea(x,y)
eff(e)

dea.plot.frontier(x,y,txt=TRUE)

# The following only works if the package FEAR is installed, it does
# not have to be loaded.
FE <- FALSE
if ( "FEAR" \%in\% .packages(TRUE) ) FE <- TRUE
# For bootstrap for real NREP should be at least 2000
# Run the following lines a couple of times and see how the bootstrap
# frontier changes from one run to the next. Try the same with
# NREP=2000
if (FE) { b <- dea.boot(x,y, NREP=100) }

# bias corrected frontier
if (FE) dea.plot.frontier(b$eff.bc*x, y, add=TRUE, lty="dashed")
# outer 95% confidence interval frontier for uncorrected frontier
if (FE) dea.plot.frontier((b$conf.int[,1]+b$bias)*x, y, add=TRUE, lty="dotted")
}


\keyword{ efficiencies }
\keyword{ bootstrap }
