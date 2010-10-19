# $Id: deaUtil.R 80 2010-10-20 08:58:43Z Lars $

efficiency <- function(object)  {
# Returnerer efficencer som et array
#   unlist(object@eff)
	e <- as.matrix(object$eff)
   if ( object$ORIENTATION == "in" )  {
      colnames(e) <- "E"
   } else if ( object$ORIENTATION == "out" )  {
      colnames(e) <- "F"
   } else if ( object$ORIENTATION == "graph" )  {
      colnames(e) <- "G"
   }
   if ( !is.null(names(object$eff)) )  {
      rownames(e) <- names(object$eff)
   }
   if ( object$TRANSPOSE == TRUE ) {
      e <- t(e)
   }
   return(e)
} ## efficiency

# effi <- efficiency
eff <- efficiency





print.Farrell  <- function(x, digits=4, ...)  {
#   a <- cbind("Efficiens"=x$eff)
    a <- x$eff
#   a <- x@eff
   print(a, digits=digits, ...)
   invisible(a)
} ## print.Farrell



summary.Farrell <- function(object, digits=4, ...)  {
   cat("Efficiency:\n")
   print.Farrell(object,digits=digits,...)
   cat("Weights (lambda):\n")
   x <- object$lambda
   xx <- round(x, digits)
   # xx <- format(unclass(x), digits=digits)
   if (any(ina <- is.na(x))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(x) < 1e-9) ) 
      xx[i0] <- sub("0", ".", xx[i0])
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(object)
   # printSpMatrix(Matrix(object$lambda),digits=4, col.names=T,...)
}  ## summary.Farrell



# returns peers, i.e. numbers for units with positive lambda,
# efficient units to be compared to
peers <- function(object, NAMES=FALSE)  {
#    if ( object$TRANSPOSE ) {
#      print("Colnames i lambda")
#      print(colnames(object$lambda))
#    } else {
#      print("Rownames i lambda")
#      print(rownames(object$lambda))
#    }
   if ( class(object) != "Farrell" && class(object) != "slack" )
      stop(paste("Object is not of class 'Farrell' (or 'slack');",
             "you might have used FAST=TRUE in 'dea'"))
   if ( object$TRANSPOSE ) {
	   lam <- object$lambda
   } else {
      lam <- t(object$lambda)
   }

   # Fjern foranstillet L_ eller L i søjlenavne for lambda
   if ( "L_" %in% substr(rownames(lam),1,2) )  {
      rownames(lam) <- substring(rownames(lam),3)
   }
   if ( "L" %in% substr(rownames(lam),1,1) )  {
      rownames(lam) <- substring(rownames(lam),2)
   }

   peer <- which(rowSums(lam,na.rm=TRUE)>0, 1:dim(lam)[1])
   # print("Firms der er peers:")
   # print(peer)

   bench <- matrix(NA,nrow=length(peer), ncol=dim(lam)[2])
   # Saet firms navne som raekkenavne
   colnames(bench) <- colnames(lam)

#   for ( i in 1:length(peer) ) {  # for hver peer
#      yper <- which(lam[peer[i],]>0)
#      bench[i,yper] <- peer[i]
#   }

   maxj = 0  # det storste antal peers for en firm
   for ( i in 1:dim(lam)[2] )  {  # for hver firm
      # Hvis firm er uden for teknologi maengde er der ingen peers: next
      if ( sum(lam[peer,i],na.rm=TRUE) == 0 ) next
      pe <- which(lam[peer,i]>0)
      bench[1:length(pe),i] <- peer[pe]
      maxj <- max(maxj,length(pe))
   }
   bench <- bench[1:maxj,,drop = FALSE]

   if (NAMES & !is.null(names(object$eff)) )  {
      bench_ <- matrix(rownames(lam)[bench], ncol=dim(bench)[2])
      # print(bench_)
      colnames(bench_) <- colnames(bench)
      bench <- bench_
   }

   if ( !object$TRANSPOSE ) {
      bench <- t(bench)
   }

#   if ( FALSE && NAMES )  {
#      # Fjern evt. foranstillet "R."
#      nv <- sub("R.","",navne.Rfirms)
#      bench <- matrix(nv[bench],nrow=dim(bench)[1])
#   }

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
   }
   if ( !object$TRANSPOSE )  lam <- t(lam)
   return(lam)
}

