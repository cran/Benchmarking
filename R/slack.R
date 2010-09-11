# $Id: slack.R 72 2010-09-11 17:06:14Z Lars $

# Calculate slack at the efficient points.

# Hvis data ikke er transponeret bliver x og y transponeret og alle
# beregninger laves transponeres. Det betyder at resultat skal
# transponeres inden afslut og retur.

slack <- function(X, Y, e, XREF=NULL, YREF=NULL, FRONT.IDX=NULL, 
                  LP=FALSE)  {
   if ( !is(e,"Farrell") )  {
       stop("In call of slack: argument 'e' must be of class 'Farrell'")
   }

   RTS <- e$RTS
   if (LP) print(paste("slack:  RTS =",RTS),quote=F)
   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }
   if ( !e$TRANSPOSE ) {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
   }

   # For at undgaa afrundingsproblemer i forbindelse med slacks skal
   # in- og ouput skaleres til at vaere i omegnen af 1
   mmm <- (rowMeans(X))
   nnn <- (rowMeans(Y))
   if ( min(mmm) < 1e-4 || max(mmm) > 1e4 || 
       min(nnn) < 1e-4 || max(nnn) > 1e4 )  {
      SKALERING <- TRUE
      X <- X / mmm
      XREF <- XREF / mmm
      Y <- Y / nnn
      YREF <- YREF / nnn
   } else
     SKALERING <- FALSE

   okr <- dim(XREF)[2]
   if (LP) cat("okr =",okr,"\n")

   if ( FALSE && length(FRONT.IDX) == 0 )  {
      # Kun units med eff==1 bruges i referenceteknologi
      FRONT.IDX <- which( abs(e$eff -1) < 1e-5 )
   }
   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in method 'slack'")
      XREF <- XREF[,FRONT.IDX, drop=FALSE]
      YREF <- YREF[,FRONT.IDX, drop=FALSE]
      # XREF <- matrix(XREF[,FRONT.IDX],nrow=dim(XREF)[1])
      # YREF <- matrix(YREF[,FRONT.IDX],nrow=dim(YREF)[1])
      if (LP) cat("FRONT.IDX =",FRONT.IDX,"\n")
   }

   m = dim(X)[1]  # number of inputs
   n = dim(Y)[1]  # number of outputs
   K = dim(X)[2]  # number of units, firms, DMUs
   Kr = dim(XREF)[2]  # number of units in reference set

   if ( is.null(e$objval) )
      eff <- e$eff
   else
      eff <- e$objval

   if ( e$ORIENTATION == "in" ) {
      E <- eff
      FF <- rep(1,K)  # Naturligt at bruge F, men F er FALSE i R
   } else if ( e$ORIENTATION == "out" )  {
      FF <- eff
      E <- rep(1,K)
   } else if ( e$ORIENTATION == "graph" )  {
      E <- eff
      FF <- 1/eff
   } else {
     stop(paste("Unknown orientation in slack:", e$ORIENTATION))
   }
   if ( m != dim(XREF)[1] )  {
      print("Number of inputs must be the same in X and XREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }
   if ( n != dim(YREF)[1] )  {
      print("Number of outputs must be the same in Y and YREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }
   if ( K !=  dim(Y)[2] )  {
      print("Number of units must be the same in X and Y",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }
   if ( Kr != dim(YREF)[2] )  {
      print("Number of units must be the same in XREF and YREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }


   if ( RTS != "crs" && RTS != "add" )  {
      rlamb <- 2
   } else {
      rlamb <- 0
   }

   # Initialiser LP objekt
   lps <- make.lp(m+n +rlamb, m+n+Kr)
   name.lp(lps, paste("DEA--slack",RTS,",",e$ORIENTATION,"orientated"))

   # saet raekker i matrix med restriktioner
   for ( h in 1:m )
       set.row(lps,h, XREF[h,], (m+n+1):(m+n+Kr) )
   for ( h in 1:n)
       set.row(lps,m+h, -YREF[h,], (m+n+1):(m+n+Kr) )
   for ( h in 1:(m+n) )
       set.mat(lps,h,h,1)
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(rep(0,m+n),rep(-1,Kr)))
      set.row(lps, m+n+2, c(rep(0,m+n),rep( 1,Kr)))
   }
   set.constr.type(lps, c(rep("=",m+n), rep(">=",rlamb)))
   set.objfn(lps, rep(1, m+n), 1:(m+n))
   lp.control(lps, sense="max")

   if ( RTS == "fdh" ) {
      set.type(lps,(m+n+1):(m+n+Kr),"binary")
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
      set.type(lps, (m+n+1):(m+n+Kr),"integer")
   }
 
   if (LP) {
      print("Right hand side not sat yet")
      print(lps)
   }

   objval <- rep(NA,K)   # vector with sums of slacks
   sx <- matrix(NA,m,K)
   sy <- matrix(NA,n,K)
   lambda <- matrix(NA,Kr,K)

   cont <- lp.control(lps)
   eps <- cont$epsilon["epsint"]

   for ( k in 1:K )  { # beregn for hver enhed

       # her er der problem med afrunding, sammen med venstrsiden er
       # det ikke sikkert at restriktionen holder naar der er
       # sket afrunding i mellemregningerner.

       rhs <- c(E[k] * X[,k], -FF[k] * Y[,k])
       set.rhs(lps, rhs, 1:(m+n))
       if (LP)  {
          print(paste("Hoejresiden for firm",k))
          print(E[k] * X[,k])
          print( -FF[k] * Y[,k])
          print(lps)
       }

      status <- solve(lps)
      if ( status != 0 ) {
	      print(paste("Error in solving for firm",k,":  Status =",status), 
           quote=F)
         objval[k] <- NA
         sol <- NA
      }  else {
         objval[k] <- get.objective(lps)
         sol <- get.variables(lps)
      }
      sx[,k] <- sol[1:m]       
      sy[,k] <- sol[(m+1):(m+n)]
      lambda[,k] <- sol[(m+n+1):(m+n+Kr)]

      if (FALSE && LP && k==1)  {
         print(paste("Obj.value =",get.objective(lps)))
         cat(paste("Solutions for",k,": "))
         print(sol)
      }

   }  # loop for each firm

   # Drop soejler i lambda der alene er nuller, dvs lambda soejler skal 
   # alene vaere reference firms

   rownames(sx) <- paste("sx",1:m,sep="")
   rownames(sy) <- paste("sy",1:n,sep="")

   # For at undgaa at regneunoejagtigheder giver ikke-nul slack bliver
   # slack sat til nul hvis de er mindre end den relative noejagitghed
   # i input og output, noejagtighed der er i beregningerne.  Det
   # betyder at sum af slacks bliver sum mens objval er sum plus
   # regneunoejagtighed som giver en forskel hvis X er meget stoerre
   # eller mindre end 1.

   sx[abs(sx) < eps ] <- 0
   sy[abs(sy) < eps ] <- 0
   if ( SKALERING )  {
      sx <- sx * mmm
      sy <- sy * nnn
   } 
   sum <- colSums(sx) + colSums(sy)

   if ( length(FRONT.IDX)>0 )  {
      rownames(lambda) <- paste("L",(1:okr)[FRONT.IDX],sep="")
   } else {
      rownames(lambda) <- paste("L",1:Kr,sep="")
   }


   if (FALSE && LP)  {
   	print("Done with loop\nColumn names for lambda:")
	   print(rownames(lambda))
	   print(lambda)
   }

   if ( !e$TRANSPOSE ) {
      sx <- t(sx)
      sy <- t(sy)
      lambda <- t(lambda)
   }

   if ( FALSE && LP ) {
      print("Faerdig med slack: lambda")
	   print(colnames(lambda))
      print(lambda)
      print("slack:")
      print(paste("laengden af objval:",length(objval)))
      print(objval==0)
      # print(!slackZero)
   }

   oe <- list(slack=sum>eps, sum=sum, objval=objval, sx=sx, sy=sy, 
        lambda=lambda, RTS=e$RTS, TRANSPOSE=e$TRANSPOSE)

   class(oe) <- "slack"

	return(oe)
} # slack






print.slack  <- function(x, digits=4, ...)  {
   sum <- x$sum
   print(sum, digits=digits, ...)
   invisible(sum)
} ## print.slack



summary.slack <- function(object, digits=4, ...)  {
   cat("Efficiency and slacks:\n")
   a <- cbind("sum"=object$sum, object$sx, object$sy, lambda(object))
   print(a,...)
   cat("Weights (lambda):\n")
   x <- object$lambda
   xx <- round(x, digits)
   # xx <- format(unclass(x), digits=digits)
   if (any(ina <- is.na(x))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(x) < 1e-9) ) 
      xx[i0] <- sub("0", ".", xx[i0])
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(xx)
}

