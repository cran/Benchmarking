# $Id: dea.direct.R 85 2010-11-06 22:18:27Z Lars $


dea.direct <- function(X,Y, DIRECT, RTS="vrs", ORIENTATION="in", 
                  XREF=NULL, YREF=NULL, FRONT.IDX=NULL, SLACK=FALSE, 
                  TRANSPOSE=FALSE)  {

   .xyref.missing <- FALSE
   if ( missing(XREF) || is.null(XREF) )  {
      .xyref.missing <- TRUE
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      .xyref.missing <- TRUE && .xyref.missing
      YREF <- Y
   }
   
   transpose <- FALSE
   if ( !TRANSPOSE )  {
      transpose <- TRUE
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }

   m <- dim(X)[1]  # number of inputs
   n <- dim(Y)[1]  # number of outputs
   K <- dim(X)[2]  # number of units, firms, DMUs
   Kr <- dim(XREF)[2] # number of units,firms in the reference technology



ee <- dea(X,Y, RTS=RTS, ORIENTATION=ORIENTATION, XREF=XREF, YREF=YREF,
         FRONT.IDX=FRONT.IDX, SLACK=SLACK, DUAL=FALSE, 
         DIRECT=DIRECT, TRANSPOSE=TRUE)  


   mmd <- switch(ORIENTATION, "in"=m, "out"=n, "in-out"=m+n) 
   ob <- matrix(ee$objval,nrow=mmd, ncol=K, byrow=T)
   if ( class(DIRECT)=="matrix" && dim(DIRECT)[2] > 1 )  {
       dir <- DIRECT
   } else {
       dir <- matrix(DIRECT,nrow=mmd, ncol=K)
   }
   if ( ORIENTATION=="in" )  {
      e <- 1 - ob*dir/X
   } else if ( ORIENTATION=="out" )  {
      e <- 1 + ob*dir/Y
   } else if ( ORIENTATION=="in-out" )  {
      e <- cbind(1 - ob[1:m,,drop=FALSE]*dir[1:m,,drop=FALSE]/X, 
        1 + ob[(m+1):(m+n),,drop=FALSE]*dir[(m+1):(m+n),,drop=FALSE]/Y)
   } else {
      warning("Illegal ORIENTATION for argument DIRECT") 
   }
   if ( class(e)=="matrix" && ( dim(e)[1]==1 || dim(e)==1 ) )
      e <- c(e) 

   if ( transpose )  {
      transpose <- FALSE
      TRANSPOSE <- FALSE
   } 


   if ( !TRANSPOSE ) {
      if ( class(e)=="matrix" )
         e <- t(e)
      ee$lambda <- t(ee$lambda)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }

   ee$eff <- e
   ee$direct <- DIRECT
   ee$TRANSPOSE <- TRANSPOSE
   return(ee)
}
