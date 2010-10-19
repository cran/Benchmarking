# $Id: dea.dual.R 78 2010-10-18 22:09:28Z Lars $

# In the calculation in the method input/output matrices X and Y are
# of the order good x firms.  Ie. X, Y etc must be transformed as
# default in R is firm x good.

# DUALIN og DUALOUT eller blot DUAL.  Her er valgt DUAL saa DUAL skal
# vaere en matrix saaledes at der en raekke i DUAL for hvert input og
# hvert output, bortset fra det foerste input og det foerste output,
# dvs m-1+n-1=m+n-2 raekker.  Der skal vaere 2 soejler, foerste soejle
# er den nedre graense, og anden soejle er den oevre graense.


dea.dual <- function(X,Y, RTS="vrs", ORIENTATION="in", 
            XREF=NULL,YREF=NULL,
            FRONT.IDX=NULL, DUAL=NULL, DIRECT=NULL,
            TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  {
   # XREF, YREF determines the technology
   # FRONT.IDX index for units that determine the technology

   rts <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( missing(RTS) ) RTS <- "vrs" 
   if (LP)  print(paste("Vaerdi af 'RTS' er ",RTS),quote=FALSE)
   if ( is.real(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=FALSE)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=FALSE)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  {
      print(paste("Unknown scale of returns:", RTS))
      print("continues asssuming RTS = \"vrs\"\n")
      RTS <- "vrs"
   } 

   orientation <- c("in-out","in","out","graph")
   if ( is.real(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      print(paste("Unknown value for ORIENTATION:",ORIENTATION),quote=F)
      ORIENTATION <- "in"
      print(paste("Continues with ORIENTATION =",ORIENTATION),quote=F)
   }

   if ( (!RTS %in% c("vrs","drs","crs","irs")) || 
              !ORIENTATION %in% c("in","out") )
      stop("dea.dual at the moment does not work for \"fdh\" or \"add\"")


   .xyref.missing <- FALSE
   if ( missing(XREF) )  {
      .xyref.missing <- TRUE
      XREF <- X
   }
   if ( missing(YREF) )  {
      .xyref.missing <- TRUE && .xyref.missing
      YREF <- Y
   }
   
   if ( !TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
   }
   orgKr <- dim(XREF)

   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in 'eff'")
      XREF <- XREF[,FRONT.IDX, drop=FALSE]
      YREF <- YREF[,FRONT.IDX, drop=FALSE]
      # XREF <- matrix(XREF[,FRONT.IDX],nrow=dim(XREF)[1])
      # YREF <- matrix(YREF[,FRONT.IDX],nrow=dim(YREF)[1])
   }

   m = dim(X)[1]  # number of inputs
   n = dim(Y)[1]  # number of outputs
   K = dim(X)[2]  # number of units, firms, DMUs
   Kr = dim(XREF)[2]  # number of units, firms, DMUs
   okr <- orgKr[2]
   if (LP) cat("m n k kr = ",m,n,K,Kr,"\n")

   if ( m != dim(XREF)[1] )  {
      print("Number of inputs must be the same in X and XREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }
   if ( n != dim(YREF)[1] )  {
      print("Number of outputs must be the same in Y and YREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   } 
   if ( K != dim(Y)[2] )  {
      print("Number of units must be the same in X and Y",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }
   if ( Kr != dim(YREF)[2] )  {
      print("Number of units must be the same in XREF and YREF",quote=F)
      return(print("Method 'eff' stops",quote=F))
   }

   if ( !is.null(DUAL) && ( !is.matrix(DUAL) ||
               dim(DUAL)[1] != (m+n-2) || dim(DUAL)[2] != 2 ) ) {
      print("DUAL must be a (m+n-2) x 2 matrix with lower and upper bounds for restrictions")
      stop("dea.dual aborts")
   }


   if ( RTS == "vrs" )  { 
      rlamb <- 2
   } else if ( RTS == "drs" || RTS == "irs" )  {
      rlamb <- 1
   } else if ( RTS == "crs" )  {
      rlamb <- 0
   } else {
      stop(paste("Unknown value for RTS in 'dea.dual':", RTS))
   }



if ( !missing(DUAL) && !is.null(DUAL) )  {
   # Make matrix for dual restrictions.
   # Foerst for input hvis der er mere end 1 input
   if ( m > 1 )  {
      Udiag <- diag(1,m-1)
      DL <- rbind(DUAL[1:(m-1),1], -Udiag)
      DU <- rbind(-DUAL[1:(m-1),2], Udiag)
      if (LP)  {
         print("DL"); print(DL)
         print("DU"); print(DU)
      }
      ADin <- cbind(DL,DU)
   } else ADin <- NULL
   if (LP)  {
      print("ADin"); print(ADin)
   }

   # og saa restriktioner for output hvis der er mere end 1 output
   if ( n > 1 )  {
      Udiag <- diag(1,n-1)
      DL <- rbind(DUAL[m:(m+n-2),1], -Udiag)
      DU <- rbind(-DUAL[m:(m+n-2),2], Udiag)
      ADout <- cbind(DL,DU)
   } else ADout <- NULL

   nulln <- matrix(0,nrow=n,ncol=2*(m-1))
   nullm <- matrix(0,nrow=m,ncol=2*(n-1))
   AD <- cbind(rbind(ADin,nulln),rbind(nullm,ADout))
   # AD <- t(AD)
   if (LP) {
      print("AD"); print(AD)
   }
} else {
  # AD <- matrix(0, nrow=m+n, ncol=2*(m+n-2))
  AD <- NULL
}

if ( is.null(AD) )  {
   restr <- 0
} else {
   restr <- 2*(m-1 +n-1)
}


   # Initialiser LP objekt
   lps <- make.lp(1+Kr+restr, m+n+rlamb )
   name.lp(lps, paste(ifelse(is.null(AD),"Dual","DualAC"),ORIENTATION,
             RTS,sep="-"))
   # if ( LP==TRUE ) print(lps)

   # saet soejler i matrix med restriktioner, saet 0'er for den foerste
   # raekke for den skal alligevel aendres for hver firm.
   for ( h in 1:m ) 
      set.column(lps,h, c( 0, -XREF[h,], AD[h,]))
   # if ( LP || !is.null(LPK) ) print(lps)
   for ( h in 1:n)
       set.column(lps,m+h, c( 0, YREF[h,], AD[m+h,]))
   # restriktioner paa lambda
   objgamma <- NULL
   if ( rlamb > 0 )  {
      set.column(lps, m+n+1, c(0,rep(-1,Kr)), 1:(Kr+1))
      objgamma <- -1
   }
   if ( rlamb > 1 )  {
      set.column(lps, m+n+2, c(0,rep( 1,Kr)), 1:(Kr+1))
      objgamma <- c(-1,1)
   }
   if ( rlamb == 1 && RTS == "irs" )  {
      set.column(lps, m+n+1, c(0,rep(1,Kr)), 1:(Kr+1))
      objgamma <- 1
   }

   set.constr.type(lps, rep("<=", 1+Kr+restr))
 
   if ( ORIENTATION == "in" )  {
      lp.control(lps, sense="max")
      set.rhs(lps, 1, 1)
   } else if ( ORIENTATION == "out" )  {
      lp.control(lps, sense="min")
      set.rhs(lps, -1, 1)
      if ( !is.null(objgamma) ) objgamma <- -objgamma
   } else
     stop("In 'dea.dual' for ORIENTATION use only 'in' or 'out'")

   if ( !is.null(CONTROL) )  {
      lp.control(lps,CONTROL)
   }
#    lp.control(lps,scaling=c("geometric","equilibrate","dynupdate"),
#                 simplextype=c("primal") )


# if ( LP || !is.null(LPK) ) print(lps)

   eff <- rep(NA,K)   # vector for the final efficiencies
   u <- matrix(NA,m,K)   # vector for the final efficiencies
   v <- matrix(NA,n,K)   # vector for the final efficiencies
   objval <- rep(NA,K)   # vector for the final efficiencies
   sol <- matrix(NA, 1 + sum(dim(lps)), K)
   gamma <- NULL
   if (rlamb > 0)  {
      gamma <- matrix(NA,rlamb,K)
   }



   for ( k in 1:K ) { # Finds the efficiencies for each unit
      # object function and first collumn
      if ( ORIENTATION == "in" )  {
         objrow <- c(rep(0,m),Y[,k], objgamma)
         if (LP) { 
            print(lps); 
            print("objgamma:") 
            print(objgamma) 
            print("objrow:")
            print(objrow)
         }
         set.objfn(lps, objrow)
         set.row(lps, 1, c(X[,k]),1:m)
      }  else  {
         objrow <- c(X[,k], rep(0,n), objgamma)
         set.objfn(lps, objrow)
         set.row(lps, 1, c(-Y[,k]),(m+1):(m+n))
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
      }  else {
         objval[k] <- get.objective(lps)
         eff[k] <- objval[k]

         losning <- get.variables(lps)
         # sol[,k] <- get.variables(lps)
         sol[,k] <- get.primal.solution(lps)
         u[,k] <- losning[1:m]
         v[,k] <- losning[(m+1):(m+n)]
         if ( rlamb > 0 )  {
            gamma[,k] <- losning[(m+n+1):(m+n+rlamb)]
         }
      }


   	if (LP && status==0) {
         print(paste("Objval, firm",k))
         print(objval[k])
         print("Solution")
         print(sol)
         print("Dual values:")
         print(u[,k])
         print(v[,k])
         print("get.variables")
         print(get.variables(lps))
         print("Primal solution")
         print(get.primal.solution(lps))
         print("Dual solution:")
         print(get.dual.solution(lps))
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         print(paste("Model",k,"(",name.lp(lps),")"))
         print(lps)
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }
   } #    for ( k in 1:K )


   # undgaa afrundingsfejl i e naar den er taet ved 1.
   eff[abs(1-eff) < 1e-5] <- 1

   rownames(u) <- paste("u",1:m,sep="")
   rownames(v) <- paste("v",1:n,sep="")

   if ( !TRANSPOSE )  {
      u <- t(u)
      v <- t(v)
   }

   oe <- list(eff=eff, objval=objval, RTS=RTS,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE,
              u=u, v=v, gamma=gamma, sol=sol)
#   class(oe) <- "Farrell"


   return (oe)
} ## dea.dual



########################

# X = x
# Y = y
# RTS="crs"
# ORIENTATION="in"
# XREF=NULL
# YREF=NULL
# FRONT.IDX=NULL
# DUAL=NULL
# TRANSPOSE=FALSE
# LP=F
# CONTROL=NULL
# LPK=NULL

