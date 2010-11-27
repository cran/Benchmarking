# $Id: deaUtil.R 97 2010-12-02 23:27:26Z Lars $


efficiencies <- function( object, ... )  {
    UseMethod( "efficiencies" )
}
eff <- function( object, ... )  {
    UseMethod( "efficiencies" )
}
eff.add <- function( object, ... )  {
    UseMethod( "efficiencies" )
}


# default method
efficiencies.default <- function( object, ... )  {
   return( object$eff )
}
eff.default <- function( object, ... )  {
   return( object$eff )
}



efficiencies.Farrell <- function(object, type="Farrell", ...)  {
# Returnerer efficencer som et array
   if ( type == "Farrell" )
     return(object$eff)
   else if ( type == "Shephard" )
     return(1/object$eff)
   else
     warning("Unknown type:", type)
   #   e <- as.matrix(object$eff)
   #   if ( object$ORIENTATION == "in" )  {
   #      colnames(e) <- "E"
   #   } else if ( object$ORIENTATION == "out" )  {
   #      colnames(e) <- "F"
   #   } else if ( object$ORIENTATION == "graph" )  {
   #      colnames(e) <- "G"
   #   }
   #   if ( !is.null(names(object$eff)) )  {
   #      rownames(e) <- names(object$eff)
   #   }
   #   if ( object$TRANSPOSE == TRUE ) {
   #      e <- t(e)
   #   }
   #   return(e)
} ## efficiency



eff.Farrell <- function(object, type="Farrell", ...)  {
      return( efficiencies.Farrell(object, type, ...) )
}





print.Farrell  <- function(x, digits=4, ...)  {
#   a <- cbind("Efficiens"=x$eff)
    a <- x$eff
#   a <- x@eff
   print(a, digits=digits, ...)
   invisible(a)
} ## print.Farrell



summary.Farrell <- function(object, digits=4, ...)  {
   eps <- 1e-6
   eff <- object$eff
   cat("Summary of efficiencies\n")
   cat("The technology is", object$RTS,"and",object$ORIENTATION,
       "orientated efficiency\n")
   cat("Number of firms with efficiency==1 are",
      sum(abs(eff-1) < eps), 
      "\nMean efficiency ", format(mean(object$eff),digit=3), "\n---" )
   if ( object$ORIENTATION!="out" )  {
      minE <- min(eff)
      minE <- floor( 10 * minE ) / 10
      dec <- seq(from=minE, to=1, by=.1)
      Estr <- "<= E <"
      Eeff <- "      E ==1   "
      n <- length(dec)
      estr <- rep(NA,n)
      for ( i in 1:(n-1) )
         estr[i] <- paste(dec[i],Estr,dec[i+1],"  ",sep="")
      estr[n-1] <- paste(estr[n-1]," ")
      estr[n] <- Eeff
      antal <- rep(NA,n)
      for ( i in 1:(n-1) )
         antal[i] <- sum(dec[i]-eps <= eff & eff < dec[i+1]-eps)
      antal[n] <- sum(abs(eff-1) < eps)
   } else {
      maxF <- max(eff)
      maxF <- ceiling( 10 * maxF ) / 10
      dec <- seq(from=1, to=maxF, by=.1)
      Estr <- "< F =<"
      Eeff <- "F ==1   "
      n <- length(dec)
      if ( n > 10 )  {
         dec_  <- c(1,1.1,1.2,1.3,1.5,2.0,5.0,10.0,100.0,Inf)
         n <- length(dec_) 
         while ( n>1 && dec_[n-1] > maxF )
            n <- n - 1 
         dec <- dec_[1:n]
      }
      estr <- rep(NA,n)
      estr[1] <- paste("    ",Eeff)
      for ( i in 2:n )
         estr[i] <- paste(format(dec[i-1],digits=2,width=3),Estr,
                          format(dec[i],digits=3,width=3),"  ",sep="")
      antal <- rep(NA,n)
      antal[1] <- sum(abs(eff-1) < eps)
      for ( i in 2:n )
         antal[i] <- sum(dec[i-1]-eps <= eff & eff < dec[i]-eps)
   }
   andel <- antal/sum(!is.na(eff))

   a <- cbind(antal , 100*andel)
   dimnames(a) <- list("  Eff range"=estr,c( "#", "%"))
   print(a,digits=c(2,3),quote=F,...)
   invisible(object)
}  ## summary.Farrell




# returns peers, i.e. numbers for units with positive lambda,
# efficient units to be compared to
peers <- function(object, NAMES=FALSE)  {
   #  if ( object$TRANSPOSE ) {
   #    print("Colnames i lambda")
   #    print(colnames(object$lambda))
   #  } else {
   #    print("Rownames i lambda")
   #    print(rownames(object$lambda))
   #  }
   if ( class(object) != "Farrell" && class(object) != "slack" )
      stop(paste("Object is not of class 'Farrell' (or 'slack');",
             "you might have used FAST=TRUE in 'dea'"))
   if ( object$TRANSPOSE ) {
	   lam <- t(object$lambda)
   } else {
      lam <- object$lambda
   }

   # Fjern foranstillet L_ eller L i søjlenavne for lambda
   if ( "L_" %in% substr(colnames(lam),1,2) )  {
      colnames(lam) <- substring(colnames(lam),3)
   }
   if ( "L" %in% substr(rownames(lam),1,1) )  {
      colnames(lam) <- substring(colnames(lam),2)
   }

   # peers har positiv vaerdi af lambda i deres soejle
   peer <- which(colSums(lam,na.rm=TRUE)>0)  #  , 1:dim(lam)[2])
   # print("Firms der er peers:")
   # print(peer)

   # Matrix til i hver raekke at holde raekkens/firmaets peers
   bench <- matrix(NA, nrow=dim(lam)[1], ncol=length(peer))
   # Saet firms navne som soejlenavne
   rownames(bench) <- rownames(lam)

   maxj = 0  # det storste antal peers for en firm
   for ( i in 1:dim(lam)[1] )  {  # for hver firm
      # Hvis firm er uden for teknologi maengden er der ingen peers: next
      if ( sum(lam[i,peer],na.rm=TRUE) == 0 ) next
      # Hvem er peers for firm i
      pe <- which(lam[i,peer]>0)
      bench[i,1:length(pe)] <- peer[pe]
      maxj <- max(maxj,length(pe))
   }
   # Der er hoejst maxj peers for en firm
   bench <- bench[,1:maxj,drop = FALSE]

   # Skal der navne i matricen bench med peers i steder for blot numre
   if (NAMES & (!is.null(rownames(lam)) || !is.null(names(object$eff))))  {
      bench_ <- matrix(rownames(lam)[bench], nrow=dim(bench)[1])
      # print(bench_)
      rownames(bench_) <- rownames(bench)
      bench <- bench_
   }

   if ( object$TRANSPOSE ) {
      bench <- t(bench)
   }

   # if ( FALSE && NAMES )  {
   #    # Fjern evt. foranstillet "R."
   #    nv <- sub("R.","",navne.Rfirms)
   #    bench <- matrix(nv[bench],nrow=dim(bench)[1])
   # }

   ## colnames(bench) <- navne.firms
   return(bench)
} ## peers



# Returns lambda-values for peers for each unit
peerslambda <- function(object)  {
   lambda <- object$lambda
   if (object$TRANSPOSE) {
	   lambda <- t(lambda)
   }
   bench <- array(NA,dim=c(2,dim(lambda)[1],dim(lambda)[2]))
   maxj = 0
   for ( h in 1:dim(lambda)[2] ) {
	j = 0
   	for ( i in 1:dim(lambda)[1] ) {
	   	if ( lambda[i,h] > 0 ) {
               j = j+1
               bench[1,j,h] = i
               bench[2,j,h] = lambda[i,h]
   		}
	   	maxj = max(maxj,j)
   	}
   }
   bench <- bench[,1:maxj,]
   return(bench)
} ## peerslambda


print.peers  <- function(x, ...)  {
   a <- peers(x)
   print(a,...)
   invisible(a)
} ## print.cost.opt



get.number.peers  <-  function(object)  {
   if ( object$TRANSPOSE ) {
	   lam <- object$lambda
   } else {
      lam <- t(object$lambda)
   }
   # Fjern foranstillet L i søjlenavne for lambda
   if ( "L" %in% substr(rownames(lam),1,1) )  {
      rownames(lam) <- substring(rownames(lam),2)
   }
   peer <- which(rowSums(lam, na.rm=TRUE)>0)
   names(peer) <- NULL
   number <- rowSums(lam[peer,]>0, na.rm=TRUE)
   cbind(peer,"#"=number)
}  # get.number.peers



get.which.peers <- function(object, N=1:length(object$eff))  {
   if ( object$TRANSPOSE ) {
	   lam <- object$lambda
   } else {
      lam <- t(object$lambda)
   }
   p <- apply(object$lambda[,N,drop=FALSE]>0,2,which)
   p0 <- p[lapply(p,length) > 0]
   return(p0)
}  # get.which.peers




lambda.print  <- function(x, KEEPREF=FALSE, ...)  {
   if ( x$TRANSPOSE ) {
      lam <- x$lambda
   } else {
      lam <- t(x$lambda)
   }
   # print(class(lam))
   if (!KEEPREF && dim(lam)[2]>1 ) {
      lam <- lam[rowSums(as.matrix(lam))>0,]
   }
   xx <- format(unclass(lam), digits=4)
   if (any(ina <- is.na(lam))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(lam) < 1e-9) ) 
      xx[i0] <- sub("0.0000", ".", xx[i0])
   if ( !x$TRANSPOSE )
      xx <- t(as.matrix(xx))
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(x)
   # printSpMatrix(Matrix(lam),digits=4, col.names=T,...)
   # invisible(lam)
} ## print.lambda



lambda <- function(object, KEEPREF=FALSE)  {
   if ( object$TRANSPOSE ) {
      lam <- object$lambda
   } else {
      lam <- t(object$lambda)
   }
   if (!KEEPREF && dim(lam)[2]>1 ) {
      lam <- lam[rowSums(lam, na.rm=TRUE)>0,,drop=FALSE]
   } else if (!KEEPREF && dim(lam)[2]==1 ) {
      lam <- lam[lam>0,,drop=FALSE]
   }

   if ( !object$TRANSPOSE )  lam <- t(lam)
   return(lam)
}



# Calculate excess input or output
excess <- function(object, X=NULL, Y=NULL)  {
   if ( class(object) != "Farrell" )
      stop("Only works for object of class/type 'Farrell'",
           "as output from dea and like functions")
   if ( is.null(object$direct) && is.null(X) && is.null(Y) )
      stop(paste("Either X or Y is needed in the arguments for",
             "objects with no direction"))

   e <- object$objval

   if ( is.null(object$direct) )  {
      # no direction, must be Farrell so direction is set by X or Y
      if ( object$ORIENTATION == "in" && !is.null(X) )
         ex <-  X * (1-e)
      else if ( object$ORIENTATION == "out" && !is.null(Y) )
         ex <-  Y * (e-1)
      else if ( object$ORIENTATION == "graph"&& !is.null(X)&& !is.null(Y) )
         ex <- cbind((1-e)*X, (1/e-1)*Y )
      else # ( is.null(dir) )
         stop("X/Y missing for ORIENTATION =", object$ORIENTATION )
   } else {
       if ( class(object$direct) == "matrix" )  {
          ex <- apply(object$direct,2,"*",e)
       } else {
           dir <- matrix(object$direct, nrow=length(e), 
                           ncol=length(object$direct), byrow=TRUE )
           ex <- e * dir
      }
   }
   # Afrund til 0 hvis ex er naer 0
   eps <- sqrt(.Machine$double.eps)
   ex[abs(ex)<eps] <- 0
   return(ex)
} # excess




eladder <- function(n, X, Y, RTS="vrs", ORIENTATION="in")  {
  idx <- NULL
  elad <- rep(NA, dim(X)[1])
  for ( i in 1:length(elad))  {
    # if (LP) print("FRONT.IDX")
    # if (LP) print(idx)
    # Brug FRONT.IDX for at kunne bruge de oprindelige indeks i X og Y
    e <- dea(X[n,,drop=F],Y[n,,drop=F], RTS=RTS, ORIENTATION=ORIENTATION,
             XREF=X,YREF=Y, FRONT.IDX=idx)
    # if (LP) print(paste("Eff =",eff(e)))
    # if (LP) print(paste("Peers =",peers(e)))
    elad[i] <- e$eff
    # Array nr. for den stoerste værdi af lambda
    p <- which(max(e$lambda)==e$lambda)
    # firm number for array number p, firm numbers follow L in the colnames
    str <- substring(colnames(e$lambda)[p],2)
    suppressWarnings(ip <- as.integer(str))
    # nok en firm der ikke laengere skal indgaa i referenceteknologien
    if ( is.na(ip) )  {
       # det er et navn/en streng saa den skal laves om til nr.
       str0 <- substring(str,2)
       navne <- rownames(X)
       if ( is.null(navne) )
           navne <- rownames(Y) 
       ip <- which( navne %in% str0 )
    }
    # saa er ip et tal
    idx <- c(idx,-ip)
  }
  return(list(eff=elad,peer=-idx))
}



eladder.plot <- function(elad, peer, TRIM=NULL)  {
   if ( !is.null(TRIM) & !is.numeric(TRIM) )
       stop("TRIM must be an integer")
   if ( is.null(TRIM) )  {
      TRIM <- 0
      for ( i in 1:length(peer) )  {
         TRIM <- max(TRIM, nchar(toString(peer[i])))
      }
   }
   linje <- ifelse(TRIM==1,2,TRIM^(1/1.3))
   opar <- 
       par(mar=c(linje+2,4.1,4.1,2.1))
   plot(elad, xaxt="n", xlab="", ylab="Efficiency")
   mtext("Most influential peers", side=1, line=linje+.5)
   if ( class(peer) == "character" || class(peer) == "factor" )  {
      axis(1, at=1:length(peer),
            labels=strtrim(peer,TRIM), las=ifelse(TRIM>1,2,0) )
   } else {
      axis(1, at=1:length(peer), labels=peer, las=ifelse(TRIM>1,2,0) )
   }
   abline(v=which(elad==1))
   abline(h=1)
   par(opar)
}

