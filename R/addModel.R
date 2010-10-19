# $Id: addModel.R 76 2010-10-12 13:45:37Z lo $

# Additive model, corresponds to eqs. 4.34-4.38 in Cooper et al., 2007 
dea.add <- function(X, Y, RTS="vrs", TRANSPOSE=FALSE,
        XREF=NULL, YREF=NULL, FRONT.IDX=NULL, LP=FALSE)  {
	if ( TRANSPOSE ) {
		K <- dim(X)[2]
	} else {
		K <- dim(X)[1]
	}
	e <- list(eff=rep(1,K), objval=NULL, RTS=RTS, 
        ORIENTATION="in", TRANSPOSE=TRANSPOSE)
   class(e) <- "Farrell"
	sl <- slack(X,Y,e, XREF=XREF, YREF=YREF, FRONT.IDX=FRONT.IDX, LP=LP)

	# Make slack the sum of slacks, not a logical variable 
	# wether there is slack or not
# 	if (TRANSPOSE)
#      sl$sum2 <- colSums(sl$sx) + colSums(sl$sy)
# 	else 
#      sl$sum2 <- rowSums(sl$sx) + rowSums(sl$sy)

	return(sl)
}  # dea.aff

