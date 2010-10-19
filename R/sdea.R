# $Id: sdea.R 75 2010-10-08 20:26:55Z Lars $


# Calculates super efficiency
sdea <- function(X,Y, RTS="vrs", ORIENTATION="in", TRANSPOSE=FALSE,
                 LP=FALSE)
{
   # Input is as for the method eff
   # Antal firmaer i data er K
   if (!TRANSPOSE)  {
      X <- t(X)
      Y <- t(Y)
   }
   K = dim(X)[2]
   if (LP)  {
      print(paste("K =",K))
      print(dim(X))
      print(dim(Y))
   }
   lambda <- matrix(NA,K,K)
   # Superefficiens skal gemmes i en K-vektor som vi til en start 
   # saetter til NA.
   supereff = rep(NA,K)
   if ( !is.null(dimnames(X)[[2]]) )  {
      names(supereff) <- dimnames(X)[[2]]
   }
   rts <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.real(RTS) )  {
      RTS_ <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTS_
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) ) {
      print(paste("Unknown value for RTS:",RTS),quote=F)
      RTS <- "vrs"
      print(paste("Continues with RTS =",RTS),quote=F)
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
   for ( i in 1:K ) {
      if (LP) print(paste("=====>> Unit",i),quote=F)
      # For hver enhed i laver vi beregningen og saetter resultat paa
      # den i'te plads i supereff
      # Den forste brug af X og Y er data for den enhed der skal 
      # beregnes efficiens for, den i'te.
      # Den anden brug er ved XREF og YREF for at angive hvilken teknologi
      # der skal bruges, de definerer teknologien.
      e <- dea(X[,i,drop=FALSE], Y[,i,drop=FALSE],
         RTS,ORIENTATION, XREF=X[,-i,drop=FALSE], YREF=Y[,-i,drop=FALSE],
         TRANSPOSE=TRUE,LP=LP)
      supereff[i] <- e$eff
      # print(dim(lambda))
      # print(dim(e$lambda))
      lambda[-i,i] <- e$lambda[,1]
   }
# print("sdea: færdig med gennemløb")

# print(colnames(X))
   if ( is.null(colnames(X)) )  {
      rownames(lambda) <- paste("L",1:K,sep="")
   } else {
       rownames(lambda) <- paste("L",colnames(X),sep="_")
   }
   colnames(lambda) <- colnames(X)
   if (!TRANSPOSE)  {
      lambda <- t(lambda)
   }

   if (LP) {
      print("sdea: Om lamda, dim og lambda")
      print(dim(lambda))
      print(rownames(lambda))
      print(colnames(lambda))
      print(lambda)
   }

   # return(supereff)
   objval <- NULL
   oe <- list(eff=supereff, lambda=lambda, objval=objval, RTS=RTS,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE,
              slack=NULL, sx=NULL, sy=NULL)
   class(oe) <- "Farrell"
   return(oe)
}  # sdea

