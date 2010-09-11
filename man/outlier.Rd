\name{outlier.ap}
\alias{outlier.ap}
\alias{outlier.ap.plot}


\title{
Detection of outliers in benchmark models
}

\description{ The functions implements the Wilson (1993) outlier
detection method using only R functions.  }

\usage{
outlier.ap(X, Y, NDEL = 3, NLEN = 25, TRANSPOSE = FALSE)

outlier.ap.plot(ratio, NLEN = 25, xlab = "r", ylab = "Log ratio", 
                ..., ylim)
}


\arguments{

  \item{X}{ Input as a firms times goods matrix, see \code{TRANSPOSE}.
}

  \item{Y}{ Output as a firms times goods matrix, see
\code{TRANSPOSE}.  }

  \item{NDEL}{ The maximum number of firms to be considered as a group
of outliers.  }

  \item{NLEN}{ The number of ratios to calculate/plot for each level
or removel, the number of rows in \code{ratio} used.  }

   \item{TRANSPOSE}{Input and output matrices are treated as goods
     times firms as defaults.  If \code{TRANSPOSE=FALSE} then marices
     are firms times goods as is standard in R for statistical models}


  \item{ratio}{ The \code{ratio} component from the list as output
from \code{outlier.ap}.} 

  \item{xlab}{ Label for the x-axis.  }

  \item{ylab}{ Label for the y-axis }

  \item{ylim}{ The y limits \code{(y1, y2)} of the plot, an
array/vector of length 2.  }

 \item{\dots}{Usual options for the methods \code{plot} and
 \code{lines}.}
 
}


\details{ Uses the function \code{det} to calculate
\eqn{R^{(i)}_{\min}}. }


\value{

\item{ratio}{}

\item{imat}{}

\item{r0}{}

}


\references{ Wilson, P.W. (1993), Detecing outliers in deterministic
nonparametric frontier models with multiple outputs, \emph{Journal of
Business and Economic Statistics} 11, 319-323.  

Wilson, P.W. (2008), "FEAR 1.0: A Software Package for Frontier
Efficiency Analysis with R," \emph{Socio-Economic Planning Sciences}
42, 247--254

}

\author{
Lars Otto \email{lo@foi.dk}
}

\note{
The function \code{outlier.ap} is extremely
slow and for \code{NDEL} larger than 3 or 4 it might be advisable to
use the function \code{ap} from the package FEAR.

The name of the returned components are the same as for \code{ap} in the package FEAR
}

\seealso{
The function ap in the package FEAR  \code{\link{help}}.
}

\examples{

n <- 25
x <- matrix(rnorm(n))
y <- .5 + 2.5*x + 2*rnorm(25)
tap <- outlier.ap(x,y,NDEL=3)
print(cbind(tap$imat,tap$rmin),na.print="",digit=2)
outlier.ap.plot(tap$ratio)

}

