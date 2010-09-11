% $Id: make.merge.Rd 72 2010-09-11 17:06:14Z Lars $
\name{make.merge}
\alias{make.merge}

\title{ Make an aggregation matrix to perform mergers of firms}

\description{ Make an aggregation matrix to perform mergers of
 firms. The matrix can be post multiplied (matrix multiplication) to
 input and output matrices to make merged input and output matrices.
 }

\usage{
make.merge(grp, nFirm = NULL, X = NULL)
}


\arguments{

  \item{grp}{ list of length Kg for Kg firms after mergers; each
component of the list is a list with the firm numbers going into this
merger.  }

  \item{nFirm}{ Number of firms before the mergers }

 \item{X}{A matrix for input or output where the rows corresponds to
the number of firms at the outset }

}


\details{ Either \code{nFirm} or \code{X} must be present; if both are
present then \code{nFirm} must be equal to the number of rows in
\code{X}, the number of firms.  }

\value{ Returns an aggregation matrix of dimension Kg times K where
rows corresponds to new merged firms and colloms are 1 for firms to be
included and 0 for firms excluded in this merged firm.  }


\note{ Whe \code{X} is an input matrix of dimension \code{K x m},
\code{K} firms and \code{m} inputs, such that if \code{M <-
make.merge(gr,K)} then \code{M \%*\% X} is the input matrix for the
merged firms. }

\author{ Peter Bogetoft and Lars Otto \email{lo@foi.dk} }

\seealso{ \code{dea.merge} }

\examples{
# To merge firms 1,3, and 5; and 2 and 4 of 7 firms into 2 new firms
# the aggregation matrix is M. Not all firms are involved in a merger.
M <- make.merge(list(c(1,3,5), c(2,4)),7)
print(M)
}

