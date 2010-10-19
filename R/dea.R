# $Id: dea.R 78 2010-10-18 22:09:28Z Lars $

# DEA beregning via brug af lp_solveAPI. Fordelene ved lp_solveAPI er
# færre kald fra R med hele matricer for hver firm og dermed skulle
# det gerne være en hurtigere metode.  Måske er det også lettere at
# gennemskue hvad der bliver gjort en gang for alle og hvad der bliver
# ændret ved beregning for hver firm.

# Option FAST=TRUE giver en meget hurtigere beregning af efficienser,
# men tilgængæld bliver der IKKE gemt de beregnede lambdaer.  Det
# betyder bl.a. at der ikke kan findes peers for de enkelte firms.


dea  <-  function(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL,YREF=NULL,
         FRONT.IDX=NULL, SLACK=FALSE, DUAL=TRUE, DIRECT=NULL,
         TRANSPOSE=FALSE, FAST=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  {
   # XREF, YREF determines the technology
   # FRONT.IDX index for units that determine the technology

   # TRANSPOSE the restriction matrix is transposed, default (TRUE) is
   # number of inputs/ouputs times number of firms, i.e. by default
   # goods are rows.

   # In the calculation in the method input/output matrices X and Y
   # are of the order good x firms.  Ie. X, Y etc must be transformed
   # as default in R is firm x good.

   if ( FAST ) { 
      DUAL=FALSE; # SLACK=FALSE; 
      # print("When  FAST then neither DUAL nor SLACK") 
   }

   rts <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.real(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  {
      print(paste("Unknown scale of returns:", RTS))
      print("continuees asssuming RTS = \"vrs\"\n")
      RTS <- "vrs"
   } 

   orientation <- c("in-out","in","out","graph")
   if ( is.real(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   if ( !(ORIENTATION %in% orientation) ) {
      print(paste("Unknown value for ORIENTATION:",ORIENTATION),quote=F)
      ORIENTATION <- "in"
      print(paste("Continues with ORIENTATION =",ORIENTATION),quote=F)
   }

   .xyref.missing <- FALSE
   if ( missing(XREF) || is.null(XREF) )  {
      .xyref.missing <- TRUE
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      .xyref.missing <- TRUE && .xyref.missing
      YREF <- Y
   }
   
   if ( !TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }
   orgKr <- dim(XREF)

   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in 'dea'")
      XREF <- XREF[,FRONT.IDX, drop=FALSE]
      YREF <- YREF[,FRONT.IDX, drop=FALSE]
      # XREF <- matrix(XREF[,FRONT.IDX],nrow=dim(XREF)[1])
      # YREF <- matrix(YREF[,FRONT.IDX],nrow=dim(YREF)[1])
   }
   rNames <- colnames(XREF)
   if ( is.null(rNames) & !is.null(colnames(YREF)) )
      rNames <- colnames(YREF)

   m = dim(X)[1]  # number of inputs
   n = dim(Y)[1]  # number of outputs
   K = dim(X)[2]  # number of units, firms, DMUs
   Kr = dim(XREF)[2] # number of units,firms in the reference technology
   oKr <- orgKr[2]
   if ( !is.null(DIRECT) )  {
      if ( class(DIRECT)=="matrix" ) {
         md <- dim(DIRECT)[1]
         Kd <- dim(DIRECT)[2]
      } else {
         md <- length(DIRECT)
         Kd <- 0
      }
   } else {
      Kd <- 0
   }
   if (LP) cat("m n K Kr = ",m,n,K,Kr,"\n")
   if (LP & !is.null(DIRECT) ) cat("md, Kd =",md,Kd,"\n") 

   if ( m != dim(XREF)[1] )  {
      print("Number of inputs must be the same in X and XREF",quote=F)
      return(print("Method 'dea' stops",quote=F))
   }
   if ( n != dim(YREF)[1] )  {
      print("Number of outputs must be the same in Y and YREF",quote=F)
      return(print("Method 'dea' stops",quote=F))
   } 
   if ( K != dim(Y)[2] )  {
      print("Number of units must be the same in X and Y",quote=F)
      return(print("Method 'dea' stops",quote=F))
   }
   if ( Kr != dim(YREF)[2] )  {
      print("Number of units must be the same in XREF and YREF",quote=F)
      return(print("Method 'dea' stops",quote=F))
   }

   if ( !is.null(DIRECT) & ORIENTATION=="graph" )
         stop("DIRECT cannot not be used with ORIENTATION=\"graph\"")
 
   if ( !is.null(DIRECT) & length(DIRECT) > 1 )  {
      if ( ORIENTATION=="in" & md!=m )
         stop("Length of DIRECT must be the number of inputs")
      else if ( ORIENTATION=="out" & md!=n )
         stop("Length of DIRECT must be the number of outputs")
      else if ( ORIENTATION=="in-out" & md!=m+n )
         stop("Length of DIRECT must be the number of inputs plus outputs")
      if ( class(DIRECT)=="matrix" & (Kd>0 & Kd!=K) )
         stop("Number of firms in DIRECT must equal firms in X and Y") 
   }
   if ( !is.null(DIRECT) & length(DIRECT) == 1 )  {
      if ( ORIENTATION=="in" & length(DIRECT)!=m )
         DIRECT <- rep(DIRECT,m)
      else if ( ORIENTATION=="out" & length(DIRECT)!=n )
         DIRECT <- rep(DIRECT,n)
      else if ( ORIENTATION=="in-out" & length(DIRECT)!=m+n )
         DIRECT <- rep(DIRECT,m+n)
   }

   if ( RTS != "crs" && RTS != "add" )  {
      rlamb <- 2
   } else 
      rlamb <- 0


   # Initialiser LP objekt
   lps <- make.lp(m+n +rlamb,1+Kr)
   name.lp(lps, paste("Dea",ORIENTATION,RTS,sep="-"))

   # saet raekker i matrix med restriktioner, saet 0'er for den foerste
   # soejle for den skal alligevel aendres for hver firm.
   for ( h in 1:m )
       set.row(lps,h, c(0,-XREF[h,]))
   for ( h in 1:n)
       set.row(lps,m+h, c(0,YREF[h,]))
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(0,rep(-1,Kr)))
      set.row(lps, m+n+2, c(0,rep( 1,Kr)))
   }

   if ( RTS == "fdh" ) {
      set.type(lps,2:(1+Kr),"binary")
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
 
   if ( !is.null(DIRECT) & Kd==0 )  {
      # print(Kd)
      # print(DIRECT)
      # Samme retning for alle enheder
      if ( ORIENTATION=="in" )
         set.column(lps, 1, c(1,-DIRECT),0:m)
      else if ( ORIENTATION=="out" )
         set.column(lps, 1, c(1,-DIRECT),c(0,(m+1):(m+n)))
      else if ( ORIENTATION=="in-out" )
         set.column(lps, 1, c(1,-DIRECT),0:(m+n))
   }

   set.objfn(lps, 1,1)
   set.constr.type(lps, rep(">=",m+n+rlamb))
   if ( ORIENTATION %in% c("in","graph") )  {
      lp.control(lps, sense="min")
   } else if ( ORIENTATION == "out" )  {
      lp.control(lps, sense="max")
   } else if ( ORIENTATION == "in-out" & !is.null(DIRECT) )  {
      lp.control(lps, sense="max")
   } else  
     stop("In 'dea' for ORIENTATION use only 'in', 'out', or 'graph'")

   if ( !is.null(DIRECT) )  {
      lp.control(lps, sense="max")
   }

   if ( !is.null(CONTROL) )  {
      lp.control(lps,CONTROL)
   }

   if ( ORIENTATION == "graph" )  {
      oe <- graphEff(lps, X, Y, XREF, YREF, RTS, FRONT.IDX, rlamb, oKr, 
                          TRANSPOSE, SLACK,FAST,LP) 
      return(oe)
   }

   objval <- rep(NA,K)   # vector for the final efficiencies
   if ( FAST ) {
     lambda <- NULL
     primal <- NULL
     dual <- NULL
   } else {
      lambda <- matrix(NA, nrow=Kr, ncol=K) # lambdas one column per unit
      rownames(lambda) <- rNames
      colnames(lambda) <- colnames(X)
      if (DUAL) {
         dual   <- matrix(NA, nrow=sum(dim(lps))+1, ncol=K) # 
         primal <- matrix(NA, nrow=sum(dim(lps))+1, ncol=K) # solutions
         colnames(dual) <- colnames(X)
         colnames(primal) <- colnames(X)
      } else {
         primal <- NULL
         dual <- NULL
      }
   }  

   # The loop for each firm
   for ( k in 1:K)  {
      # Af en eller anden grund saetter set.column ogsaa vaerdi for
      # kriteriefunktion og hvis der ikke er nogen vaerdi bliver den
      # automatisk sat til 0.  Derfor maa 1-tallet for
      # kriteriefunktionen med for denne soejle og det er raekke 0.

      if ( is.null(DIRECT) )  {
         if ( ORIENTATION == "in" )  {
            set.column(lps, 1, c(1,X[,k]),0:m)
            set.rhs(lps, Y[,k], (m+1):(m+n))
         } else {
            set.column(lps, 1, c(1,-Y[,k]),c(0,(m+1):(m+n)))
            set.rhs(lps, -X[,k], 1:m)
         }
      } else {
         # print(Kd)
         # print(DIRECT)
         set.rhs(lps, c(-X[,k],Y[,k]), 1:(m+n))
         if ( Kd > 1 )  {
         # retning for enheden
            if ( ORIENTATION=="in" )
               set.column(lps, 1, c(1,-DIRECT[,k]),0:m)
            else if ( ORIENTATION=="out" )
               set.column(lps, 1, -DIRECT[,k],0:n)
            else if ( ORIENTATION=="in-out" )
               set.column(lps, 1, c(1,-DIRECT[,k]),0:(m+n))
         }
      }
      if ( LP )  print(paste("Firm",k), quote=FALSE)
      if ( LP && k == 1 )  print(lps)
      status <- solve(lps)
      if ( status != 0 ) {
        if (status == 2) {
	        print(paste("Firm",k,"not in the technology set"), quote=F)
        } else {
	        print(paste("Error in solving for firm",k,":  Status =",status), 
             quote=F)
        }
        objval[k] <- NA
        sol <- NA
      }  else {
         objval[k] <- get.objective(lps)
         if ( !FAST ) sol <- get.variables(lps)
      }
      if ( !FAST )  {
         lambda[,k] <- sol[2:(1+Kr)]
         if ( DUAL )  {
            primal[,k] <- get.primal.solution(lps)
            dual[,k] <- get.dual.solution(lps)
         }
      }

   	if (LP && status==0) {
         print(paste("Objval, firm",k))
         print(get.objective(lps))
         print("Solution/varaibles")
         print(get.variables(lps))
         print("Primal solution")
         print(get.primal.solution(lps))
         print("Dual solution:")
         print(get.dual.solution(lps))
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }
   }  # loop for each firm

   # Afrund efficiencer der naermest er 1 til 1 saa der ikke er
   # afrundingsfejl; noget stoerre end sqrt(.Machine$double.eps)
   e <- objval
   lpcontr <- lp.control(lps)
   eps <- lpcontr$epsilon["epsint"]
   e[abs(e-1) < eps] <- 1
   if ( !is.null(dimnames(X)[[2]]) )  {
      names(e) <- dimnames(X)[[2]]
   }

#   if ( ORIENTATION == "in" )  {
#      names(e) <- "E"
#   } else if ( ORIENTATION == "out" )  {
#      names(e) <- "F"
#   } else if ( ORIENTATION == "graph" )  {
#      names(e) <- "G"
#   }

   if ( FAST ) { 
      return(e)
      stop("Her skulle vi ikke kunne komme i 'dea'")
   }
   if (LP) print("Forbi retur fra FAST")

   if ( is.null(rownames(lambda)) )  {
      if ( length(FRONT.IDX)>0 )  {
         rownames(lambda) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
      } else {
         rownames(lambda) <- paste("L",1:Kr,sep="")
      }
   } else {
       rownames(lambda) <- paste("L",rownames(lambda),sep="_")
   }

   sign <- NULL
   if ( ORIENTATION == "out" ) sign <- -1 else sign <- 1

   if ( DUAL )  {
     ux <- sign*dual[2:(1+m),,drop=FALSE] 
     vy <- sign*dual[(2+m):(1+m+n),,drop=FALSE] 
     rownames(ux) <- paste("u",1:m,sep="")
     rownames(vy) <- paste("v",1:n,sep="")
     if ( rlamb > 0 ) 
        gamma <- dual[(1+m+n+1):(1+m+n+rlamb),,drop=FALSE]
     else
        gamma <- NULL
     
   } else {
     ux <- vy <- NULL
   }
   if (LP) print("DUAL faerdig")  
   
   if ( !TRANSPOSE ) {
      lambda <- t(lambda)
      if (DUAL)  {
         ux <- t(ux)
         vy <- t(vy)
         primal <- t(primal)
         dual <- t(dual)
         if ( !is.null(gamma) ) gamma <- t(gamma)
      }
   }

   oe <- list(eff=e, lambda=lambda, objval=objval, RTS=RTS,
              primal=primal, dual=dual, ux=ux, vy=vy, gamma=gamma,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE
              # ,slack=slack_, sx=sx, sy=sy
              )
   class(oe) <- "Farrell"


   if ( SLACK ) {
      if ( !TRANSPOSE )  { # Transponer tilbage hvis de blev transponeret
         X <- t(X)
         Y <- t(Y)
         if (.xyref.missing) {
            XREF <- NULL
            YREF <- NULL
         } else {
            XREF <- t(XREF)
            YREF <- t(YREF)
         }
      }
      sl <- slack(X, Y, oe, XREF, YREF, FRONT.IDX, LP=LP)
      oe$slack <- sl$slack
      oe$sx <- sl$sx
      oe$sy <- sl$sy
      oe$lambda <- sl$lambda
      if (LP)  {
         print("slack fra slack:")
         print(sl$slack)
         print("slack efter slack:")
         print(oe$slack)
      }
   }


   return(oe)

}

