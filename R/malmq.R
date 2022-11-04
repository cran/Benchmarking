# $Id: malmq.R 247 2022-09-15 22:38:27Z X052717 $

# Beregner Malmquist indeks og dekomponering af samme ud fra to perioder

# Det er ikke forudsat at der er samme antal units i hver periode, men
# indekset bliver kun beregnet for foreningsmaengden af units

# Metoden kan bruges alene, men er ogsaa taenkt som kaldt fra en anden
# metode for haandtering af flere perioder.

 malmq <- function(X0, Y0, ID0=NULL,  X1,Y1, ID1=NULL, 
          RTS="vrs", ORIENTATION="in", SAMEREF=FALSE,
          SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
          TRANSPOSE=FALSE, FAST=TRUE, LP=FALSE, CONTROL=NULL, LPK=NULL)  
{

   # Det er underforstaaet at X'er og Y'er er matricer 'units x var'

   rts <- c("fdh","vrs","drs","crs","irs","irs2","add","fdh+","fdh++","fdh0")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.numeric(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  stop("Unknown scale of returns: ", RTS)

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop("Unknown value for ORIENTATION: ", ORIENTATION)
   }

   m0<- dim(X0)[2]  # number of inputs
   n0<- dim(Y0)[2]  # number of outputs
   K0 <- dim(X0)[1]  # number of units, firms, DMUs
   m1 <- dim(X1)[2]  # number of inputs
   n1 <- dim(Y1)[2]  # number of outputs
   K1 <- dim(X1)[1]  # number of units, firms, DMUs

   # Hvis ID'er mangler bliver de sat til 1,...,K
   if ( is.null(ID0) )  ID0 <- seq(1,K0)
   if ( is.null(ID1) )  ID1 <- seq(1,K1)

   # Laengden af ID skal svare til antal units i X og Y
   if ( K0 != length(ID0) ||  K0 != dim(Y0)[1] )
      stop("Number of units in X0 and Y0 must correspont to length of ID0")
   if ( K1 != length(ID1) ||  K1 != dim(Y1)[1] )
      stop("Number of units in X1 and Y1 must correspont to length of ID1")

    # Antal input skal vaere samme i X0 og X1
    if (m0 != m1) stop("Number of inputs must be the same in X0 and X1")
    # Antal output skal vaere samme i Y0 og Y1
    if (n0 != n1) stop("Number of outputs must be the same in Y0 and Y1")

   # Find faellesmaengden af units og tilhoerende indeks
   idlab <- intersect(ID0, ID1)
   id0 <- ID0 %in% idlab
   id1 <- ID1 %in% idlab

    # Er der gengangere, og er der samme antal i de to perioder
    if ( sum(id0) != length(unique(ID0[id0])) || sum(id1) != length(unique(ID1[id1])) || 
            sum(id0)!=sum(id1) )
        stop("Units in ID are not unique for each period")

    # Input og output for faellesmaengden af units. 
    # De samme units der skal beregnes for ellers giver det ingen mening
    # at beregne Malmquist indeks
    x0 <- X0[id0,,drop=FALSE]   
    y0 <- Y0[id0,,drop=FALSE]
    x1 <- X1[id1,,drop=FALSE]   
    y1 <- Y1[id1,,drop=FALSE]
    
    # Reference teknologi kan godt vaere bestemt af forskelige units i 
    # forskellige perioder, men det kan ogsaa vaere de samme units der 
    # bestemmer teknologien i hver periode.
    if (SAMEREF)  {
        X0 <- X0[id0,,drop=FALSE]   
        Y0 <- Y0[id0,,drop=FALSE]
        X1 <- X1[id1,,drop=FALSE]   
        Y1 <- Y1[id1,,drop=FALSE]
    }

   # Skal teknologien i en periode bestemmes af alle dem der i
   # perioden uafhaengigt af om de er i den anden periode?  Som det er
   # gjort nu er teknologien bestemt af dem der er i begge perioder.
   # Ved at bruge X0, X1 mm. som reference bestemmes teknologien af alle
   # units i paagaeldende periode, selv om der kun beregnes indeks for units
   # der er i begge perioder.

   e00_ <- dea(x0, y0, RTS=RTS, ORIENTATION=ORIENTATION, XREF=X0,YREF=Y0,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e10_ <- dea(x1, y1, RTS=RTS, ORIENTATION=ORIENTATION, XREF=X0,YREF=Y0,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e11_ <- dea(x1, y1, RTS=RTS, ORIENTATION=ORIENTATION, XREF=X1,YREF=Y1,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e01_ <- dea(x0, y0, RTS=RTS, ORIENTATION=ORIENTATION, XREF=X1,YREF=Y1,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)

    # Lav raekkefoelgen saa den svarer til den i idlab
    e00 <- rep(NA, length(idlab))
    e01 <- rep(NA, length(idlab))
    e10 <- rep(NA, length(idlab))
    e11 <- rep(NA, length(idlab))
    for( i in idlab)  {
        e00[idlab==i] <- e00_[ID0[id0]==i]
        e01[idlab==i] <- e01_[ID0[id0]==i]
        e10[idlab==i] <- e10_[ID1[id1]==i]
        e11[idlab==i] <- e11_[ID1[id1]==i]
    }
   tc <- sqrt(e10/e11 * e00/e01)  # teknisk aendring; flytning af frontier
   ec <- e11/e00                  # aendring i effektivitet
   m  <- tc * ec
   mq <- sqrt(e10/e00 * e11/e01)  ## == m;  Malmquist indeks for produktivitet

   return( list(m=m, tc=tc, ec=ec, mq=mq, id=idlab, id0=id0, id1=id1, 
                e00=e00, e10=e10, e11=e11, e01=e01) )

}
