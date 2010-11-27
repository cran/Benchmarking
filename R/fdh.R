# $Id: fdh.R 96 2010-11-27 22:47:25Z Lars $

# FDH efficiency uden brug af LP.
# Der er ingen kontrol af argumenter, den taenkes at blive kaldt fra dea.R
# der udfoerer alle kontroller


fdh <- function(X,Y, ORIENTATION="in", XREF=NULL, YREF=NULL, 
                FRONT.IDX=NULL, DIRECT=NULL, TRANSPOSE=FALSE)  {

   if ( ORIENTATION=="graph" )
      stop("ORIENTATION==\"graph\" does not work for fdh,",
            "use dea( ...,RTS=\"fdh\", ... ) ")

   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }

   if ( TRANSPOSE )  {
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
   }
   rNames <- colnames(XREF)
   if ( is.null(rNames) & !is.null(colnames(YREF)) )
      rNames <- colnames(YREF)

   K <- dim(X)[1]
   m <- dim(X)[2]
   n <- dim(Y)[2]
   Kr <- dim(XREF)[1]
   eff <- rep(NA, K)
   peer <- rep(NA, K)


# Find dominating reference firms for each of the Kr reference firms
# Aendres saa der kun dropppes firms der ikke kan give loesning,
# mindre input kan vaere relevant ved input super efficiens og mere
# output ved output super efficiens.

# Do <- list(NA,Kr)
# for ( i in 1:Kr )  {
#    # brug af rowSums er hurtigere end brug af apply eller all()
#   Do[[i]] <- which(
#    rowSums(matrix(XREF[i,,drop=FALSE],nrow=Kr,ncol=m,byrow=TRUE) >= 
#            XREF) == m
#    & 
#    rowSums(matrix(YREF[i,,drop=FALSE],nrow=Kr,ncol=n,byrow=TRUE) <= 
#            YREF) == n
#    )
#    # Hvis der er andre firms der dominerer den firmen selv fra listen 
#    # fordi den kan have slack og derfor optræde med efficens 1 og så være 
#    # peer, men det kan give probmlem med flere peers.
#    if ( length(Do[[i]]) > 1 ) {
#       i0 <- which( Do[[i]]==i )
#       Do[[i]] <- Do[[i]][-i0]
#    }
# }
# DDo <- unique(unlist(Do))


# Directional efficiency
if ( !is.null(DIRECT) )  {

   if ( class(DIRECT)=="matrix" && dim(DIRECT)[1] > 1 ) {
      if ( ORIENTATION=="in" )  {
         dirX <- DIRECT  # matrix(DIRECT, nrow=K, ncol=m)
         # dirY <- matrix(.Machine$double.xmin, nrow=K, ncol=n)
         dirY <- matrix(NA, nrow=K, ncol=n)
      } else if ( ORIENTATION=="out" )  {
         # dirX <- matrix(.Machine$double.xmin, nrow=K, ncol=m)
         dirX <- matrix(NA, nrow=K, ncol=m)
         dirY <- DIRECT  # matrix(DIRECT, nrow=K, ncol=n)
      } else if ( ORIENTATION=="in-out" )  {
         dirX <- DIRECT[,1:m,drop=FALSE]   
                 # matrix(DIRECT[,1:m], nrow=K, ncol=m)
         dirY <- DIRECT[,(m+1):(m+n), drop=FALSE]  
                 # matrix(DIRECT[,(m+1):(m+n)], nrow=K, ncol=n)
      }
   } else {
      # Her er DIRECT en vektor og derfor ens for alle firms, dvs.
      # alle raekker skal vaere ens
      if ( ORIENTATION=="in" )  {
         dirX <- matrix(DIRECT, nrow=K, ncol=m, byrow=T)
         # dirY <- matrix(.Machine$double.xmin, nrow=K, ncol=n)
         dirY <- matrix(NA, nrow=K, ncol=n)
      } else if ( ORIENTATION=="out" )  {
         # dirX <- matrix(.Machine$double.xmin, nrow=K, ncol=m)
         dirX <- matrix(NA, nrow=K, ncol=m)
         dirY <- matrix(DIRECT, nrow=K, ncol=n, byrow=T)
      } else if ( ORIENTATION=="in-out" )  {
         dirX <- matrix(DIRECT[1:m], nrow=K, ncol=m, byrow=T)
         dirY <- matrix(DIRECT[(m+1):(m+n)], nrow=K, ncol=n, byrow=T)
      }
   }

   for ( k in 1:K )  {
      # For each firm find max(XREF/X) over inputs (rows)
      # Se kun for de firmaer der dominerer firm k i den retning der 
      # ikke aendres
      xk <-  NULL  # matrix(X[k,,drop=FALSE], nrow=Kr, ncol=m, byrow=TRUE)
      yk <-  NULL  # matrix(Y[k,,drop=FALSE], nrow=Kr, ncol=n, byrow=TRUE)
      if ( ORIENTATION=="in" )  {
         yk <- matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
         idx <- rowSums(yk <= YREF) == n
      } else if ( ORIENTATION=="out" )  {
         xk <-  matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
         idx <- rowSums(xk >= XREF) == m
      } else if ( ORIENTATION=="in-out" )  {
         xk <-  matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
         yk <-  matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
         idx <- rowSums(xk >= XREF) == m & rowSums(yk <= YREF) == n
      }
      if ( is.null(xk) ) 
         xk <-  matrix(X[k,,drop=FALSE], nrow=sum(idx), ncol=m, byrow=TRUE)
      else
         xk <-  xk[idx,,drop=FALSE]

      if ( is.null(yk) )  
         yk <-  matrix(Y[k,,drop=FALSE], nrow=sum(idx), ncol=n, byrow=TRUE)
      else 
         yk <-  yk[idx,,drop=FALSE]

      allDir <- cbind( (xk-XREF[idx,])/dirX[rep(k,sum(idx)),], 
                       (YREF[idx,]-yk)/dirY[rep(k,sum(idx)),])
      minDir <- apply(allDir,1,min, na.rm=TRUE)
      eff[k] <- max(minDir)
      # der er kun gjort plads til een peer per firm
      peer[k] <- (1:Kr)[idx][which.max(minDir)]
   }
   # we only need lambda to be able to call peers() to get peers.
   lam <- matrix(0, nrow=K, ncol=Kr)
   for (k in 1:K)  {
       lam[k, peer[k]] <- 1
   }
   e <- list(eff=eff, objval=eff, peers=peer, lambda=lam, RTS="fdh",
             direct=DIRECT, ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
   class(e) <- "Farrell"
   return(e)
}  # if !is.null(DIRECT)




# Find efficiency when compared to each dominating firm; the actual 
# efficiency is then then min/max over the dominating firms.
if ( ORIENTATION=="in" )  {
   for ( k in 1:K )  {
      # For each firm find max(XREF/X) over inputs (rows)
      yk <- matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
      idx <- rowSums(yk <= YREF) == n
      if ( sum(idx) == 0 )  {
         # Der er ingen loesning, eff er NA
         eff[k] <- Inf
         peer[k] <- NA
         next
      } 
      maxIn <- apply( XREF[idx,,drop=FALSE] / 
        matrix(X[k,], nrow=sum(idx), ncol=m, byrow=TRUE)
                   , 1, max )
      eff[k] <- min(maxIn)
      # der er kun gjort plads til een peer
      peer[k] <- (1:Kr)[idx][which.min(maxIn)]
   }
} else {
   for ( k in 1:K )  {
      xk <- matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
      idx <- rowSums(xk >= XREF) == m
      if ( sum(idx) == 0 )  {
         # Der er ingen loesning, eff er NA
         eff[k] <- -Inf
         peer[k] <- NA
         next
      } 
      minOut <- apply( YREF[idx,,drop=FALSE] / 
        matrix(Y[k,,drop=FALSE], nrow=sum(idx), ncol=n, byrow=TRUE)
                      , 1, min)
      eff[k] <- max(minOut)
      peer[k] <- (1:Kr)[idx][which.max(minOut)]
   }
}

# we only need lambda to be able to call peers() to get peers.
lam <- matrix(0, nrow=K, ncol=Kr)
for (k in 1:K)  {
    lam[k, peer[k]] <- 1
}
# print(lam)


e <- list(eff=eff, objval=eff, peers=peer, lambda=lam, RTS="fdh", 
          ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
class(e) <- "Farrell"

return(e)
}  # fdh function

