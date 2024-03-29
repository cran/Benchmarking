# $Id: dea.R 257 2023-06-07 09:51:25Z larso $

# DEA beregning via brug af lp_solveAPI. Fordelene ved lp_solveAPI er
# faerre kald fra R med hele matricer for hver unit og dermed skulle
# det gerne vaere en hurtigere metode.  Maaske er det ogsaa lettere at
# gennemskue hvad der bliver gjort en gang for alle og hvad der bliver
# aendret ved beregning for hver unit.

# Option FAST=TRUE giver en meget hurtigere beregning af efficienser,
# men tilgaengaeld bliver der IKKE gemt de beregnede lambdaer.  Det
# betyder bl.a. at der ikke kan findes peers for de enkelte units.


dea  <-  function(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL,YREF=NULL,
         FRONT.IDX=NULL, SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
         TRANSPOSE=FALSE, FAST=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  
{
    # XREF, YREF determines the technology
    # FRONT.IDX index for units that determine the technology

    # In the calculation in the method input/output matrices X and Y
    # are of the order good x units.  

    # TRANSPOSE the restriction matrix is transposed, For TRUE then X and Y are
    # matrices of dimension inputs/ouputs times number of units, i.e.
    # goods are rows, and therefore X, Y etc must be transformed
    # as default in R is unit x good.

    if ( FAST ) { 
        DUAL=FALSE; # SLACK=FALSE; 
        # print("When  FAST then neither DUAL nor SLACK") 
    }

    #----------
    # Tjek af RTS og ORIENTATION
    rts <- c("fdh","vrs","drs","crs","irs","irs2","add","fdh+","fdh++","fdh0", "vrs+")
    if ( missing(RTS) ) RTS <- "vrs" 
    if ( is.numeric(RTS) )  {
        if (LP) print(paste("Number '",RTS,"'",sep=""),quote=F)
        RTStemp <- rts[1+RTS] # the first, fdh, is number 0
        RTS <- RTStemp
        if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=F)
    }
    RTS <- tolower(RTS)
    if ( !(RTS %in% rts) )  stop("Unknown scale of returns:", RTS)

    orientation <- c("in-out","in","out","graph")
    if ( is.numeric(ORIENTATION) )  {
        ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
        ORIENTATION <- ORIENTATION_
    }
    ORIENTATION <- tolower(ORIENTATION)
    if ( !(ORIENTATION %in% orientation) ) {
        stop("Unknown value for ORIENTATION:",ORIENTATION)
    }

    #----------
    # Tjek af input
    # Hvis data er en data.frame saa tjek om det er numerisk data og lav
    # dem i saa fald om til en matrix
    X <- tjek_data(X)
    Y <- tjek_data(Y)

    .xyref.missing <- FALSE
    if ( missing(XREF) || is.null(XREF) )  {
        .xyref.missing <- TRUE
        XREF <- X
    }
    if ( missing(YREF) || is.null(YREF) )  {
        .xyref.missing <- TRUE && .xyref.missing
        YREF <- Y
    }
   
    XREF <- tjek_data(XREF)
    YREF <- tjek_data(YREF)

    if ( TRANSPOSE )  {
        X <- t(X)
        Y <- t(Y)
        XREF <- t(XREF)
        YREF <- t(YREF)
        if ( !is.null(DIRECT) & is(DIRECT, "matrix") )
            DIRECT <- t(DIRECT)
    }

    #------------
    
    # Tjek af dimensioner
    orgKr <- dim(XREF)
    if ( length(FRONT.IDX) > 0 )  {
        if (LP) print("FRONT.IDX")
        if (LP) print(FRONT.IDX)
        if ( !is.vector(FRONT.IDX)  & !(ifelse(is.matrix(FRONT.IDX), dim(FRONT.IDX)[2]==1, FALSE) ))
            stop("FRONT.IDX is not a vector or collumn matrix in 'dea'")
        XREF <- XREF[FRONT.IDX,, drop=FALSE]
        YREF <- YREF[FRONT.IDX,, drop=FALSE]
    }
    rNames <- rownames(XREF)
    if ( is.null(rNames) & !is.null(rownames(YREF)) )
        rNames <- rownames(YREF)
    if ( is.null(rNames) )
        rNames <- 1:dim(XREF)[1]

    m <- dim(X)[2]  # number of inputs
    n <- dim(Y)[2]  # number of outputs
    K <- dim(X)[1]  # number of units, units, DMUs
    Ky <- dim(Y)[1]  
    Kr <- dim(XREF)[1] # number of units,units in the reference technology
    oKr <- orgKr[1]    # number of units in reference before use of FRONT.IDX
    if ( !is.null(DIRECT) )  {
        if ( is(DIRECT, "matrix") ) {
            md <- dim(DIRECT)[2]
            Kd <- dim(DIRECT)[1]
        } else {
            md <- length(DIRECT)
            Kd <- 0
        }
    } else {
        Kd <- 0
    }
    if (LP) cat("m n K Kr = ",m,n,K,Kr,"\n")
    if (LP & !is.null(DIRECT) ) cat("md, Kd =",md,Kd,"\n") 

    if (m!=dim(XREF)[2])  stop("Number of inputs must be the same in X and XREF")
    if (n!=dim(YREF)[2])  stop("Number of outputs must be the same in Y and YREF")
    if (K!=Ky)  stop("Number of units must be the same in X and Y")
    if (Kr!=dim(YREF)[1])  stop("Number of units must be the same in XREF and YREF")

    if ( !is.null(DIRECT) & all(DIRECT=="min") & ORIENTATION=="graph" )
        # Kaldet kommer fra 'mea' og stop vil vise 'dea' kaldet som mea laver; derfor call.=FALSE
        stop("The option 'ORIENTATION=\"graph\"' cannot be used for 'mea'", call.=FALSE)
 
    if ( !is.null(DIRECT) & length(DIRECT) > 1 )  {
      if ( ORIENTATION=="in" & md!=m )
         stop("Length of DIRECT must be the number of inputs")
      else if ( ORIENTATION=="out" & md!=n )
         stop("Length of DIRECT must be the number of outputs")
      else if ( ORIENTATION=="in-out" & md!=m+n )
         stop("Length of DIRECT must be the number of inputs plus outputs")
      if ( is(DIRECT, "matrix") & (Kd>0 & Kd!=K) )
         stop("Number of units in DIRECT must equal units in X and Y") 
    }
    if ( !is.null(DIRECT) & length(DIRECT) == 1 )  {
      if ( ORIENTATION=="in" & length(DIRECT)!=m )
         DIRECT <- rep(DIRECT,m)
      else if ( ORIENTATION=="out" & length(DIRECT)!=n )
         DIRECT <- rep(DIRECT,n)
      else if ( ORIENTATION=="in-out" & length(DIRECT)!=m+n )
         DIRECT <- rep(DIRECT,m+n)
    }

    #------------
    
    if ( RTS=="fdh" && ORIENTATION!="graph" && !FAST && DUAL==FALSE
            && all(DIRECT!="min") )  {
        e <- fdh(X,Y, ORIENTATION=ORIENTATION, XREF=XREF, YREF=YREF, 
               FRONT.IDX=FRONT.IDX, DIRECT=DIRECT, TRANSPOSE=FALSE, oKr)
        if ( SLACK )  {
            warning("Run 'slack(X, Y, e)' to get slacks")
        }
        return(e)
    }
    
    #------------

    if ( RTS=="fdh++" )  {
        e <- dea.fdhPlus(X, Y, ORIENTATION=ORIENTATION,
            XREF=XREF, YREF=YREF, FRONT.IDX=FRONT.IDX, DIRECT=DIRECT, 
            param=param, TRANSPOSE=FALSE, oKr)
        if ( SLACK )  {
            warning("Run 'slack(X, Y, e)' to get slacks")
        }
        return(e)
    }
    
    #------------

    if ( RTS != "crs" && RTS != "add" )  {
        rlamb <- 2
    } else {
        rlamb <- 0
    }

    #------------

    # Lav LP objekt og saet alt der ikke aendres ved ny firm.
    # Initialiser LP objekt
    lps <- make.lp(m+n +rlamb,1+Kr)
    # Saet lp options
    lp.control(lps,
        scaling=c("range", "equilibrate", "integers")  # default scalering er 'geometric'
    )                   # og den giver ikke altid tilfredsstillende resultat;
                        # curtisreid virker i mange tilfaelde slet ikke
    if ( is.null(DIRECT) ) dirStreng<-"" else dirStreng<-"Dir"
    name.lp(lps, paste(ifelse(is.null(DIRECT)|(DIRECT!="min"), "Dea", "Mea"),
                            ORIENTATION,RTS,dirStreng,sep="-"))
    
    if (!missing(CONTROL)) set_control(lps, CONTROL)

    # saet raekker i matrix med restriktioner, saet 0'er for den foerste
    # soejle for den skal alligevel aendres for hver unit.
    # Bemaerk X og Y transponeres implicit naar de saettes i lp_solveAPI; 
    # soejler i X og Y saettes som raekker i lp_solveAPI.
    # Foerste 'm' raekker med input
    for ( h in 1:m )
       set.row(lps,h, c(0,-XREF[,h]))
    # Foelgende 'n' raekker med output
    for ( h in 1:n)
       set.row(lps,m+h, c(0,YREF[,h]))
    # restriktioner paa lambda
    if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(0,rep(-1,Kr)))
      set.row(lps, m+n+2, c(0,rep( 1,Kr)))
    }

    # Saet restriktioner for lambda, dvs. for RTS
    if ( RTS == "fdh" || RTS == "fdh0" ) {
      set.type(lps,2:(1+Kr),"binary")
      set.rhs(lps,-1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
    } else if ( RTS == "vrs" )  {
      set.rhs(lps, c(-1,1), (m+n+1):(m+n+2))
    } else if ( RTS == "vrs+" )  {
    # param: (lav, hoej, sum lav, sum hoej)
      # Saet parametrene low og high
      if ( is.null(param) )  {
         param <- c(.5, 2.)
      }
      if ( length(param) == 1 )  {
         low <- param
         high <- 1+(1-param)
      } else {
         low <- param[1]
         high <- param[2]
      }
      if (length(param) == 4)  {
            # print("Graenser for sum af lambda sat")
            set.rhs(lps, c(-param[4], param[3]), (m+n+1):(m+n+2))
        } else {
            set.rhs(lps, c(-1,1), (m+n+1):(m+n+2))
        }
        param <- c(low=low, high=high)
        set.semicont(lps, 2:(1+Kr))
        set.bounds(lps, lower=rep(low,Kr), upper=rep(high,Kr), columns=2:(1+Kr))
    } else if ( RTS == "drs" )  {
        set.rhs(lps, -1, m+n+1)
        delete.constraint(lps, m+n+2)
        rlamb <- rlamb -1
    } else if ( RTS == "irs" )  {
        set.rhs(lps, 1, m+n+2)
        delete.constraint(lps, m+n+1)
        rlamb <- rlamb -1
    } else if ( RTS == "irs2" )  {
        set.rhs(lps, 1, m+n+2)
        delete.constraint(lps, m+n+1)
        rlamb <- rlamb -1
        set.semicont(lps, 2:(1+Kr))
        set.bounds(lps, lower=rep(1,Kr), columns=2:(1+Kr))
    } else if ( RTS == "add" )  {
        set.type(lps,2:(1+Kr),"integer")
    } else if ( RTS == "fdh+" )  {
        # Saet parametrene low og high
        if ( is.null(param) )  {
           param <- .15
        }
        if ( length(param) == 1 )  {
           low <- 1-param
           high <- 1+param
        } else {
           low <- param[1]
           high <- param[2]
        }
        param <- c(low=low, high=high)
        set.rhs(lps, c(-high, low), (m+n+1):(m+n+2))
        add.SOS(lps,"lambda", 1,1, 2:(1+Kr), rep(1, Kr))
    }

    if ( !is.null(DIRECT) & Kd<=1 & all(DIRECT!="min") )  {
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

    if ( !is.null(DIRECT) )  {
      # Ved super efficiency for directional kan loesning vaere
      # negativ saa loensingen skal kunne vaere negativ, dvs. objval
      # kan vaere negativ.
      set.bounds(lps,lower=-Inf, columns=1)
    }

    set.objfn(lps, 1,1)
    # Baade in- og output modeller skal formuleres med ">="
    set.constr.type(lps, rep(">=",m+n+rlamb))
    if ( ORIENTATION %in% c("in","graph") )  {
      lp.control(lps, sense="min")
    } else if ( ORIENTATION == "out" )  {
       lp.control(lps, sense="max")
    } else if ( ORIENTATION == "in-out" & !is.null(DIRECT) )  {
       lp.control(lps, sense="max")
    } else  
      stop("In 'dea' for ORIENTATION use only 'in', 'out', 'graph', or 'in-out (only for DIRECT)")

    # Ved brug af directional efficiency er der altid tale om et max-problem
    if ( !is.null(DIRECT) )  {
      lp.control(lps, sense="max")
    }

    if ( ORIENTATION == "graph" )  {
        oe <- graphEff(lps, X, Y, XREF, YREF, RTS, FRONT.IDX, rlamb, oKr, 
                          param=param, TRANSPOSE, SLACK, FAST, LP) 
        # delete.lp(lps)
        return(oe)
    }

    objval <- rep(NA,K)   # vector for the final efficiencies
    if ( FAST ) {
        lambda <- NULL
        primal <- NULL
        dual <- NULL
    } else {
        lambda <- matrix(NA, nrow=K, ncol=Kr) # lambdas one column per unit
        rownames(lambda) <- rownames(X)
        colnames(lambda) <- rNames
        if (DUAL) {
            dual   <- matrix(NA, nrow=K, ncol=sum(dim(lps))+1) # 
            primal <- matrix(NA, nrow=K, ncol=sum(dim(lps))+1) # solutions
            rownames(dual) <- rownames(X)
            rownames(primal) <- rownames(X)
        } else {
            primal <- NULL
            dual <- NULL
        }
    }

    if ( !is.null(DIRECT) & all(DIRECT == "min") )  {
        directMin <- TRUE
        if ( ORIENTATION=="in" )  {
            directMatrix <- matrix(NA, nrow=K, ncol=m)
        } else if ( ORIENTATION=="out" )  {
            directMatrix <- matrix(NA, nrow=K, ncol=n)
        } else if ( ORIENTATION=="in-out" )  {
            directMatrix <- matrix(NA, nrow=K, ncol=m+n)
        }
    } else {
        directMin <- FALSE
    }

    # The loop for each unit
    for ( k in 1:K)  {
      if ( LP )  print(paste("Unit",k," -------------------"), quote=FALSE)
 
      # Af en eller anden grund saetter set.column ogsaa vaerdi for
      # kriteriefunktion og hvis der ikke er nogen vaerdi bliver den
      # automatisk sat til 0.  Derfor maa 1-tallet for
      # kriteriefunktionen med for denne soejle og det er raekke 0.


      if ( directMin )  {
         # Saet hoejreside for enhedens input og output
         set.rhs(lps, c(-X[k,],Y[k,]), 1:(m+n))

         # Find retningen og saet foerste soejle til den
         if ( ORIENTATION=="in" )  {
            DIRECT <- minDirection(lps, m, n, ORIENTATION, LP=LP)
            set.column(lps, 1, c(1,-DIRECT),0:m)
         } else if ( ORIENTATION=="out" )  {
            DIRECT <- minDirection(lps, m, n, ORIENTATION, LP=LP)
            set.column(lps, 1, c(1,-DIRECT), c(0,(m+1):(m+n)) )
         } else if ( ORIENTATION=="in-out" )  {
            DIRECT <- minDirection(lps, m, n, ORIENTATION)
            set.column(lps, 1, c(1,-DIRECT),0:(m+n))
         }
         directMatrix[k,] <- DIRECT
         if (LP) { print("Min DIRECT:"); print(DIRECT) }

         # Check om DIRECT er 0, hvis den er nul gaa til naeste unit
         lpcontr <- lp.control(lps)
         eps <- sqrt(lpcontr$epsilon["epsint"])
         if (LP) print(paste("eps for minDirectin",eps), quote=FALSE)
# Daarligt test for om direction er 0, tager ikke hensyn at X'er og
# Y'er indbyrdes kan vaere forskellig stoerrelse
         if ( ORIENTATION=="in" )  {
            deltaDir <- DIRECT/( X[k,] + .Machine$double.xmin )
         } else if ( ORIENTATION=="out" )  {
            deltaDir <- DIRECT/( Y[k,] + .Machine$double.xmin )
         } else if ( ORIENTATION=="in-out" )  {
            deltaDir <- DIRECT/( c(X[k,],Y[k,]) + .Machine$double.xmin )
         }

         if ( max(DIRECT) < eps && max(abs(deltaDir)) < eps )
         {
            if (LP) print(paste("Direction 0 for unit",k))
            objval[k] <- 0
            if ( !FAST )  {
               lambda[k,] <- rep(0,Kr)
               lambda[k,k] <- 1
            }
            next  # ingen direction at gaa, tag naeste unit
         }
      }  # if ( directMin )


      if ( is.null(DIRECT) )  {
         if ( ORIENTATION == "in" )  {
            set.column(lps, 1, c(1,X[k,]),0:m)
            set.rhs(lps, Y[k,], (m+1):(m+n))
         } else {
            set.column(lps, 1, c(1,-Y[k,]),c(0,(m+1):(m+n)))
            set.rhs(lps, -X[k,], 1:m)
         }
      } else {
         # print(Kd)
         # print(DIRECT)
         set.rhs(lps, c(-X[k,],Y[k,]), 1:(m+n))
         if ( Kd > 1 )  {
         # retning for enheden
            if ( ORIENTATION=="in" )
               set.column(lps, 1, c(1,-DIRECT[k,]),0:m)
            else if ( ORIENTATION=="out" )
               set.column(lps, 1, c(1,-DIRECT[k,]), c(0,(m+1):(m+n)) )
            else if ( ORIENTATION=="in-out" )
               set.column(lps, 1, c(1,-DIRECT[k,]),0:(m+n))
         }
      }
      if ( LP && k <= 10 )  print(lps)
      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if ( status == 5 ) {
         # Numerical failure, reset basis og proev igen evt. med en 
         # anden skalering; desvaerre virker det helt forkerte resultater
         # at aendre skalering, dvs. at bruge 'dynupdate'.
         set.basis(lps, default=TRUE)
         status <- solve(lps)
      }
      if (LP)  print(paste("Status =",status))
      if ( status != 0 )  {
        if ( status == 2 || status == 3 ) {
          # At unit ikke er i teknologimaengden svarer til 'Inf', det
          # er derfor unoedvendigt at give en advarsel.
          #cat("Unit ", k, " is not in the technology set. Status = ", status,
          #    "\n", sep="")
          objval[k] <- ifelse(ORIENTATION=="in", Inf, -Inf)
        } else {
           cat("Error in solving for unit ", k, ":  Status =",status, "\n", sep="")
           objval[k] <- NA
        }
        sol <- NA
      }  else {
         objval[k] <- get.objective(lps)
         if ( !FAST ) sol <- get.variables(lps)
      }
      if ( !FAST )  {
         lambda[k,] <- sol[2:(1+Kr)]
         if ( DUAL )  {
            primal[k,] <- get.primal.solution(lps)
            dual[k,] <- get.dual.solution(lps)
         }
      }

      if (LP && status==0) {
         print(paste("Objval, unit",k))
         print(get.objective(lps))
         print("Solution/variables")
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
    }  # loop for each unit

    e <- objval

    lpcontr <- lp.control(lps)
    eps <- lpcontr$epsilon["epsint"]
    e[abs(e-1) < eps] <- 1      # 'e' taet ved 1 skal vaere 1
    if ( !is.null(dimnames(X)[[1]]) )  {
       names(e) <- dimnames(X)[[1]]
    }
    lambda[abs(lambda-1) < eps] <- 1   # taet ved 1
    lambda[abs(lambda) < eps] <- 0     # taet ved 0

    # Faerdig med at bruge lps
    # delete.lp(lps)

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
         colnames(lambda) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
      } else {
         colnames(lambda) <- paste("L",1:Kr,sep="")
      }
    } else {
       colnames(lambda) <- paste("L",rNames,sep="_")
    }

    if ( DUAL )  {
        if ( ORIENTATION == "out" ) sign <- -1 else sign <- 1
        ux <- sign*dual[,2:(1+m),drop=FALSE] 
        vy <- sign*dual[,(2+m):(1+m+n),drop=FALSE] 
        colnames(ux) <- paste("u",1:m,sep="")
        colnames(vy) <- paste("v",1:n,sep="")
        if ( rlamb > 0 ) {
            gamma <- dual[,(1+m+n+1):(1+m+n+rlamb),drop=FALSE]
		   if (rlamb==1)  {
			   colnames(gamma) <- "gamma"
			} else {
				colnames(gamma) <- paste("gamma",1:rlamb,sep="")
			}
        } else
            gamma <- NULL
     
    } else {
        ux <- vy <- NULL
    }  ##  if (DUAL)
    if (LP) print("DUAL faerdig")

    if ( directMin )  {
        DIRECT <- directMatrix
    }
 
   
    if ( TRANSPOSE ) {
        if ( is(e, "matrix") )
            e <- t(e)
        lambda <- t(lambda)
        if (DUAL)  {
            ux <- t(ux)
            vy <- t(vy)
            primal <- t(primal)
            dual <- t(dual)
            if ( !is.null(gamma) ) gamma <- t(gamma)
        }
        if ( !is.null(DIRECT) & is(DIRECT, "matrix") )
            DIRECT <- t(DIRECT)
    }  ## if (TRANSPOSE)

    oe <- list(eff=e, lambda=lambda, objval=objval, RTS=RTS,
              primal=primal, dual=dual, ux=ux, vy=vy, gamma=gamma,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE,
              # slack=NULL, sx=NULL, sy=NULL, sum=NULL, 
              param=param 
              )

    if (!is.null(DIRECT))  {
        oe$direct <- DIRECT
    }
   
    class(oe) <- "Farrell"


    if ( SLACK )  {
         if ( TRANSPOSE )  { # Transponer tilbage hvis de blev transponeret
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
        oe$sum <- sl$sum
        oe$sx <- sl$sx
        oe$sy <- sl$sy
        oe$lambda <- sl$lambda
        if (LP)  {
            print("slack fra slack:")
            print(sl$slack)
            print("slack efter slack:")
            print(oe$slack)
        }
    }  ## if (SLACK)

   return(oe)

}  # dea



# Kontrol af om data er numerisk
data_kontrol <- function(X)  {
    if ( is.null(X) )  return(TRUE)
    if (!is.data.frame(X) && !is.numeric(X))  return(FALSE)
    if ( is.data.frame(X))  {
       nc <- dim(X)[2]
       for ( i in 1:nc )  {
          if ( !is.numeric(X[,i]) )  {
             return(FALSE)
             break
          }
       }
    }
    return(TRUE)
}  # data_kontrol



# Tjek om data er matrix og hvis ikke lav om til matrix
tjek_data <- function(X)  {
    if (is.matrix(X) & is.numeric(X)) return(X)
    if (is.null(X) || !data_kontrol(X)) {
        stop("'", deparse(substitute(X)), 
             "' is not a numeric matrix (or numeric data.frame)", call.=FALSE)
    }
    if (!is.data.frame(X) || !is.numeric(X))  {
        X <- as.matrix(X)
    }
    return(X)
}  # tjek_data



# Saet option til lpsolveAPI; isaer relevant for skalering
set_control <- function(lp, CONTROL=NULL)
{
    if ( !is.null(CONTROL) )  {
        if( !is.list(CONTROL)) {
            stop( "argument 'CONTROL' must be a 'list' object")
        }
        do.call( lp.control, c( list( lprec = lp ), CONTROL ) )
    }
    return(NULL)
}  # set_control
