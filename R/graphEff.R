# $Id: graphEff.R 72 2010-09-11 17:06:14Z Lars $

# Funktion til beregning af graf efficiens.  Beregning sker via
# bisection hvor der itereres mellem mulige og ikke-mulige løsninger
# et LP problem hvor venstreside er som in- og output orienteret
# efficiens. blot er første søjle erstattet af rene nuller, og G*X og
# (1/G)*Y optræder på højresiden. Minimering af 0 så der blot søges om
# der er en mulig løsning. Da søjlen for efficiens er bar 0'er vil
# justering af efficiens ud over G ikke ske, dvs. det er kun lambdaer
# der tilpasses for at se om der er en mulig løsning.

graphEff <- function(lps, X, Y, XREF, YREF, RTS, FRONT.IDX, rlamb, oKr, 
                          TRANSPOSE, SLACK,FAST,LP) 
{
   m = dim(X)[1]  # number of inputs
   n = dim(Y)[1]  # number of outputs
   K = dim(X)[2]  # number of units, firms, DMUs
   Kr = dim(YREF)[2]  # number of units, firms, DMUs

   objval <- rep(NA,K)   # vector for the final efficiencies
   if ( FAST ) {
     lambda <- NULL
   } else {
      lambda <- matrix(NA, nrow=Kr, ncol=K) # lambdas one column per unit
   }
   set.column(lps, 1, rep(0,dim(lps)[1]))
   tol <- 1e-5
   for ( k in 1:K)  {
      if ( LP )  print(paste("Firm",k), quote=FALSE)
      # Lav bisection
      a <- 0
      b <- 2  # medfører start med G=1
      nIter <- 0
      while ( b-a > tol && nIter < 50 )  {
         G <- (a+b)/2
         set.rhs(lps, c(-G*X[,k],Y[,k]/G), 1:(m+n))
         # if ( k==1 ) print(lps)
         status <- solve(lps)
         if (LP) print(paste("G = ",G,"(",k,"); status =",status))
         if ( status == 0 ) {
            # løsning findes
            b <- G
         } else {
            a <- G
         }
         nIter <- nIter + 1
      }
      if ( status != 0 )  {
         # Hvis den sidste værdi af G ikke var mulig bruger vi den
         # øvre grænse. Det er nødvendigt med en mulig løsning for at
         # kunne få lambdaer og duale værdier.
         G <- b
	      set.rhs(lps, c(-G*X[,k],Y[,k]/G), 1:(m+n))
         status <- solve(lps)
	   }
      if (LP)  {
         print(paste("G = ",G,"(",k,"); status =",status))
         print(rlamb)
         print("Solution")
         print(get.variables(lps))
         print(lps)
      }
      objval[k] <- G
      if ( LP && k == 1 )  print(lps)
      if ( !FAST ) 
      if ( !FAST )  {
         sol <- get.variables(lps)
         lambda[,k] <- sol[2:(1+Kr)]
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
   }  # loop for each firm

   e <- objval
   # lpcontr <- lp.control(lps)
   # eps <- lpcontr$epsilon["epsint"]
   e[abs(e-1) < tol] <- 1

   if ( FAST ) { 
      return(e)
      stop("Her skulle vi ikke kunne komme i 'dea'")
   }

   if ( length(FRONT.IDX)>0 )  {
      rownames(lambda) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
   } else {
      rownames(lambda) <- paste("L",1:Kr,sep="")
   }

   primal <- dual <- NULL
   ux <- vy <- NULL

   if ( !TRANSPOSE ) {
      lambda <- t(lambda)
   }

   oe <- list(eff=e, lambda=lambda, objval=objval, RTS=RTS,
              primal=primal, dual=dual, ux=ux, vy=vy, gamma=gamma,
              ORIENTATION="graph", TRANSPOSE=TRANSPOSE
              # ,slack=slack_, sx=sx, sy=sy
              )
   class(oe) <- "Farrell"

   if ( SLACK ) {
      if ( !TRANSPOSE )  { # Transponer tilbage hvis de blev transponeret
         X <- t(X)
         Y <- t(Y)
         XREF <- t(XREF)
         YREF <- t(YREF)
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

