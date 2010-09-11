% $Id: sfa.Rd 72 2010-09-11 17:06:14Z Lars $
\name{sfa}
\alias{sfa}
\alias{lambda.sfa}
\alias{logLik.sfa}
\alias{print.sfa}
\alias{residuals.sfa}
\alias{sigma2.sfa}
\alias{sigma2u.sfa}
\alias{sigma2v.sfa}
\alias{summary.sfa}
\alias{coef.sfa}
\alias{te.sfa}
\alias{te1.sfa}
\alias{te2.sfa}
\alias{teJ.sfa}
\alias{teBC.sfa}
\alias{teMode.sfa}
\alias{te.add.sfa}


\title{Stochastic frontier estimation}

\description{
Estimate a stochastic frontier production functin using a maximum 
likelihood method.
}

\usage{
sfa(x, y, beta0 = NULL, lambda0 = 1, resfun = ebeta, 
    TRANSPOSE = FALSE, DEBUG=FALSE,..., control, hessian=2)
}


\arguments{

  \item{x}{input as a k x m matrix of observations on m inputs of
k firms; (unit x input); MUST be a matrix. No constant for the
intercept should be included in x as it is added by default.}

  \item{y}{output; k times 1 matrix (one output)  }

  \item{beta0}{optional initial parameter values}

%  \item{sigma0}{optional initial estimate of total variance}

  \item{lambda0}{optional initial ratio of variances}

  \item{resfun}{function to calculate the residuals, default is a
  linear model with an intercept. Must be called as
  \code{resfun(x,y,parm)} where \code{parm=c(beta,lambda)} or
  \code{parm=c(beta)}, and return the residuals as an array of length
  corresponding to the length of output \code{y}.}

%   \item{t.value}{Calculates the t-values and write the results with
%   the parameters; also calculates standard errors of the parameters}
% 
%   \item{std.err}{Calculates the standard errors of the parameter
%   estimates and write the results along the parameters}
% 
%   \item{vcov}{Calculates the varians-covarians matrix of the parameters}

  \item{TRANSPOSE}{Data is transposed, ie. input is now m x k matrix}

  \item{DEBUG}{Set to TRUE to get various debuggen informatin written
  on the console}

  \item{\dots}{optional parameters for the \code{ucminf} method}

  \item{control}{list of control parameters to \code{ucminf}}

  \item{hessian}{how the Hessian is delivered, see the ucminf
  documentation}
}


\details{
  The optimization is done by the R method \code{ucminf}. The efficiency terms are asumed to be half--normal distributed.
}

\value{

The values returned is the same as for \code{ucminf}, i.e. a list with components plus some espicifically relevant for sfa: 

\item{par}{The best set of parameters found \code{c(beta,lambda)}.}

\item{value}{The value of minus loglikelihood function corresponding to 'par'.}

\item{beta}{The parameters for the function}

\item{sigma2}{The estimate of the total variance}

\item{lambda}{The estimate of lambda}

\item{N}{The number of observations}

\item{df}{The  degrees of freedom for the model}

\item{residuals}{The residuals as a n times 1 matrix/vector, 
can also be obtained by \code{residuals(sfa-object)}}

\item{fitted.values}{Fitted values}

\item{vcov}{The varians-covarians matrix for all estimaed parameters
 incl. lambda}

\item{convergence}{An integer code. '0' indicates successful
          convergence.  Some of the error codes taken from
          \code{ucminf} are

'1' Stopped by small gradient (grtol). 

'2' Stopped by small step (xtol). 

'3' Stopped by function evaluation limit (maxeval). 

'4' Stopped by zero step from line search 

More codes are found in \code{\link{ucminf}}
}

\item{message}{A character string giving any additional information
    returned by the optimizer, or 'NULL'.}


\item{o}{The object returned by \code{ucminf}, for further information
on this see \code{ucminf}}

}

\references{See the method \code{ucminf} for the possible optimization
methods and further options}

\author{Lars Otto \email{lo@foi.dk}}

\note{ Calculation of technical efficiencys for each unit can be done
by the method te.sfa as shown in the examples.

\code{te.sfa(sfaObject), teBC.sfa(sfaObject)}:  Battese and Coelli (1988, 392)

\code{teMode.sfa(sfaObject), te1.sfa(sfaObject)}:  Jondrow et al. (1982, 235), conditional mode

\code{teJ.sfa(sfaObject), te2.sfa(sfaObject)}:  Jondrow et al. (1982, 235), conditional mean

\code{te.add.sfa(sfaObject)} efficiency in the additive model

The variance for the distribution of efficiencies can be calculated by
\code{sigma2u.sfa(sfaObject)}, the variance for the random errors by
\code{sigma2v.sfa(sfaObject)}, and the total variance by
\code{sigma2.sfa}.

}


\seealso{\code{\link{ucminf}}}

\examples{
# Example from the book by Coelli et al.
# d <- read.csv("c:/0work/rpack/front41Data.csv", header = TRUE, sep = ",")
# x <- cbind(log(d$capital), log(d$labour))
# y <- matrix(log(d$output))

n <- 50
x1 <- 1:50 + rnorm(50,0,10)
x2 <- 100 + rnorm(50,0,10)
x <- cbind(x1,x2)
y <- 0.5 + 1.5*x1 + 2*x2 + rnorm(n,0,1) - pmax(0,rnorm(n,0,1))
sfa(x,y)
summary(sfa(x,y))


# Estimate efficiency for each unit
o <- sfa(x,y)
te <- te.sfa(o)
te1 <- te1.sfa(o)
te2 <- te2.sfa(o)
cbind(te,te1,te2)[1:10,]

# Detailed calculations with out the method te.sfa
e <- o$residuals    # Residuals
residuals(o)[1:10]        # or

sigma2.sfa(o)       # Estimated varians
lambda.sfa(o)       # Estimated lambda


}
\keyword{sfa}
\keyword{efficiency}%
