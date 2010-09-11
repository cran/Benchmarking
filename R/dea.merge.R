# $Id: dea.merge.R 72 2010-09-11 17:06:14Z Lars $


dea.merge <- function(X, Y, M, RTS = "vrs", ORIENTATION = "in", 
    XREF = NULL, YREF = NULL,
    FRONT.IDX = NULL, TRANSPOSE = FALSE, LP = FALSE, CONTROL=NULL,
    LPK=NULL, ...)
{

# Mergers should be measured against the originale technology set,
# therefore we must use XREF og YREF otherwise it will be Xmerger and
# Ymerger that determines the technology.

.xyref.missing <- FALSE
if ( missing(XREF) || is.null(XREF) )  {
   .xyref.missing <- TRUE
   XREF <- X
}
if ( missing(YREF) || is.null(YREF) )  {
   .xyref.missing <- TRUE && .xyref.missing
   YREF <- Y
}

Xmerger <- M %*% X
Ymerger <- M %*% Y

# Potential gains
E <- dea(Xmerger, Ymerger, RTS=RTS, ORIENTATION=ORIENTATION, 
    XREF = XREF, YREF = YREF,
    FRONT.IDX = FRONT.IDX, TRANSPOSE = TRANSPOSE, FAST = TRUE, 
    LP = LP, CONTROL=CONTROL, LPK=LPK, ...)


# Individual efficiencies
e <- dea(X, Y, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE,  
      LP = LP, CONTROL=CONTROL, LPK=LPK)


# Inputs of individual firms projected on efficient frontier and
# inputs of merged firms after elimination of individual inefficiency
if ( ORIENTATION == "in" )  {
   Xeff <- diag(e) %*% X
   Xmerger_proj <- M %*% Xeff
   Ymerger_proj <- Ymerger
} else if ( ORIENTATION == "out" )  {
   Yeff <- diag(e) %*% Y
   Xmerger_proj <- Xmerger
   Ymerger_proj <- M %*% Yeff
} else if ( ORIENTATION == "graph" )  {
   Xeff <- diag(e) %*% X
   Yeff <- diag(e) %*% Y
   Xmerger_proj <- M %*% Xeff
   Ymerger_proj <- M %*% Yeff
} else  {
   stop(paste("Unknown ORIENTATION:",ORIENTATION,
              " The function dea.merge stops"))
}

# Pure gains from mergers
Estar <- dea(Xmerger_proj, Ymerger_proj, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE,  
      LP = LP, CONTROL=CONTROL, LPK=LPK)


# Learning effect
LE <- E/Estar 

# Inputs and outputs for merged firms in harmony calculation
Xharm <- diag(1/rowSums(M)) %*% Xmerger_proj
Yharm <- diag(1/rowSums(M)) %*% Ymerger_proj

# Harmony effect
HA <- dea(Xharm,Yharm, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE,  
      LP = LP, CONTROL=CONTROL, LPK=LPK) 


# Size effect
SI <- Estar/HA

em <- list(Eff=E, Estar=Estar, learning=LE, harmony=HA, size=SI)

return(em)

} 
