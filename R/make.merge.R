# $Id: make.merge.R 72 2010-09-11 17:06:14Z Lars $


make.merge <- function(grp, nFirm=NULL, X=NULL)  {
   # Opstiller aggregeringsmatrix for at danne grupperne grp ud fra X.
   # Hvad der skal merges skal angives som indeks i en liste af arrays
   # hvor hvert array er indeks for de enheder der skal indg� i en given
   # gruppe
   g <- length(grp)
   # print(g)
   if ( is.null(X) & is.null(nFirm) ) {
      stop("Either X or nFirm must be in the call to merge.matrix")
   }
   Kx <- -1
   if ( !is.null(X) ) {
      K <- Kx <- dim(X)[1]
   }
   if ( !is.null(nFirm) )
      K <- nFirm
   if ( !is.null(nFirm) & !is.null(X) & Kx != K )
      stop("nFirm must be the number of rows in X")
   Mer <- matrix(0, nrow=g, ncol=K)
   for ( i in 1:g )  {  # S�t 1-taller s�jler for dem der skal merges
	   # print(paste("Gruppe", i)) 
	   # print(grp[[i]])
     	Mer[i,grp[[i]]] = 1 
   }
   return(Mer)    # returnerer merge matrix
   # X %*% Mer    # returnerer merge input/output data
}