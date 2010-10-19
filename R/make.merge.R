# $Id: make.merge.R 81 2010-10-20 15:26:25Z Lars $


make.merge <- function(grp, nFirm=NULL, X=NULL)  {
   # Opstiller aggregeringsmatrix for at danne grupperne grp ud fra X.
   # Hvad der skal merges skal angives som indeks i en liste af arrays
   # hvor hvert array er indeks for de enheder der skal indgå i en given
   # gruppe
   if ( class(grp) == "factor" )  {
      g <- nlevels(grp)
      K <- Kg <- length(grp)
   } else {
     g <- length(grp)
     Kg <- -1
   }
   # print(g)
   if ( Kg == -1 & is.null(X) & is.null(nFirm) ) {
      stop("Either X or nFirm must be in the call to merge.matrix or grp must be a factor")
   }
   Kx <- -1
   if ( !is.null(X) ) {
      K <- Kx <- dim(X)[1]
   }
   if ( !is.null(nFirm) )
      K <- nFirm
   if ( !is.null(nFirm) & !is.null(X) & Kx != K )
      stop("nFirm must be the number of rows in X")
   if (Kg!=-1 & !is.null(nFirm) & Kg!=K )
      stop("nFirm must be the length of the facotr grp")
   if (Kg!=-1 & !is.null(X) & Kg!=Kx )
      stop("The length of the factor grp must be the number of rows in X")
   Mer <- matrix(0, nrow=g, ncol=K)
   if ( class(grp) == "factor" )  {
     for ( i in 1:g )  {  # Sæt 1-taller søjler for dem der skal merges
     	    Mer[i,as.numeric(grp)==i] <- 1 
     }
   } else {
     for ( i in 1:g )  {  # Sæt 1-taller søjler for dem der skal merges
	       # print(paste("Gruppe", i)) 
	       # print(grp[[i]])
     	    Mer[i,grp[[i]]] <- 1 
     }
   }
   return(Mer)    # returnerer merge matrix
   # X %*% Mer    # returnerer merge input/output data
}

