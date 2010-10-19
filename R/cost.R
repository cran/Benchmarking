# $Id: cost.R 75 2010-10-08 20:26:55Z Lars $

# Function to calculate minimum cost input.

# Calculatins are done with trasnposed matrices compared to R
# standards, but according to LP practice

cost.opt <- function(XREF, YREF, W, YOBS=NULL, RTS="vrs", TRANSPOSE=FALSE, 
                     LP=FALSE, LPK=NULL)  {
   if ( missing(YOBS) )  {
   	YOBS <- YREF
   }
   if (!TRANSPOSE) {
      XREF <- t(XREF)
      YREF <- t(YREF)
      W <- t(W)
      YOBS <- t(YOBS)
   }
   m = dim(XREF)[1]  # number of inputs
   n = dim(YREF)[1]  # number of outputs
   K = dim(YOBS)[2]  # number of units, firms, DMUs
   Kr = dim(XREF)[2]  # number of units, firms, DMUs

   if ( dim(W)[2] > 1 && dim(W)[2] != K )
      stop("Dimensions for W and YOBS are different")
   if ( Kr != dim(YREF)[2] )
      stop("Number of firms in XREF and YREF differ")
   if ( m != dim(W)[1] )
      stop("Number of inputs in W and XREF differ")
   if ( n != dim(YOBS)[1] )
      stop("Number of outputs in YREF and YOBS differ")

   rts <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.real(RTS) )  {
      if (LP) cat(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) cat(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   if ( !(RTS %in% rts) )  {
      print(paste("Unknown scale of returns:", RTS))
      print("continuees asssuming RTS = \"vrs\"\n")
      RTS <- "vrs"
   } 

   if ( RTS != "crs" && RTS != "add" )  {
      rlamb <- 2
   } else 
      rlamb <- 0

   lps <- make.lp(m+n +rlamb,m+Kr)
   name.lp(lps, paste("DEA cost,",RTS,"technology"))

   # saet raekker i matrix med restriktioner, saet 0'er for den foerste
   # soejle for den skal alligevel aendres for hver firm.
   dia <- diag(1,nrow=m)
   for ( h in 1:m )
       set.row(lps,h, c(dia[h,], -XREF[h,]))
   for ( h in 1:n)
       set.row(lps,m+h, c(rep(0,m), YREF[h,]))
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(rep(0,m),rep(-1,Kr)))
      set.row(lps, m+n+2, c(rep(0,m),rep( 1,Kr)))
   }

   if ( RTS == "fdh" ) {
      set.type(lps,2:(1+K),"binary")
      set.rhs(lps,-1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
   } else if ( RTS == "vrs" )  {
      set.rhs(lps, c(-1,1), (m+n+1):(m+n+2))
   } else if ( RTS == "drs" )  {
      set.rhs(lps, -1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
   } else if ( RTS == "irs" )  {
      set.rhs(lps, 1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
   } else if ( RTS == "add" )  {
      set.type(lps,2:(1+Kr),"integer")
   }

   set.objfn(lps, c(W[,1],rep(0,K)))
   set.constr.type(lps, rep(">=",m+n+rlamb))
   lp.control(lps, sense="min")

   xopt <- matrix(NA,m,K)
   lambda <- matrix(NA,nrow=Kr,ncol=K)
   cost <- matrix(NA,1,K)

   for ( k in 1:K ) {
      if ( dim(W)[2] != 1 && k > 1 ) { 
         set.objfn(lps, c(W[,k],rep(0,K)))
	   }
      set.rhs(lps, YOBS[,k], (m+1):(m+n))

      if (LP) print(lps)

      status <- solve(lps)
      if ( status != 0 ) {
	      print(paste("Error in solving for firm",k,":  Status =",status), 
           quote=F)
      }  else {
         cost[k] <- get.objective(lps)
         sol <- get.variables(lps)
         xopt[,k] <- sol[1:m]
         lambda[,k] <- sol[(m+1):(m+Kr)]
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }

   }  # for ( k in 1:K )

   rownames(lambda) <- paste("L",1:Kr,sep="")
   rownames(cost) <- rownames(YOBS)

   if (!TRANSPOSE) {
      xopt <- t(xopt)
      lambda <- t(lambda)
   }

   svar <- list("xopt"=xopt, "cost"=cost, "lambda"=lambda, 
                rts=RTS, TRANSPOSE=TRANSPOSE)
   class(svar) <- "cost.opt"
   return (svar)
} # cost.opt



print.cost.opt  <- function(x, ...)  {
   a <- cbind("Optimalt input"=x$x)
   print(a,...)
   invisible(a)
} ## print.cost.opt


summary.cost.opt <- function(object, ...)  {
   cat("Optimal input:\n")
   print.cost.opt(object)
   cat("Cost:\n")
   print(object$cost,...)
   cat("Weights (lambda):\n")
   x <- object$lambda
   xx <- format(unclass(x), digits=4)
   if (any(ina <- is.na(x))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(x) < 1e-9) ) 
      xx[i0] <- sub("0.0000", ".", xx[i0])
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(object)
   # printSpMatrix(Matrix(object$lambda),digits=4, col.names=T,...)
   # print(object$lambda,digits=4)
}  ## summary.cost.opt


