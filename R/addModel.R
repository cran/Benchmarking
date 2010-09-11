# $Id: addModel.R 72 2010-09-11 17:06:14Z Lars $

# Additive model, corresponds to eqs. 4.34-4.38 in Cooper et al., 2007 
add.dea <- function(X, Y, RTS="vrs", TRANSPOSE=FALSE,...)  {
	if ( TRANSPOSE ) {
		K <- dim(X)[2]
	} else {
		K <- dim(X)[1]
	}
	e <- list(eff=rep(1,K), objval=NULL, RTS=RTS, 
        ORIENTATION="in", TRANSPOSE=TRANSPOSE)
   class(e) <- "Farrell"
	sl <- slack(X,Y,e,...)

	# Make slack the sum of slacks, not a logical variable 
	# wether there is slack or not
# 	if (TRANSPOSE)
#      sl$sum2 <- colSums(sl$sx) + colSums(sl$sy)
# 	else 
#      sl$sum2 <- rowSums(sl$sx) + rowSums(sl$sy)

	return(sl)
}  # aff.dea
