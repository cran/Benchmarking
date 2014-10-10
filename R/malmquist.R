# $Id$

# Beregner Malmquist indeks for enhederne ID over tidspunkterne i TIME

# Det forudsættes at der ikke er huller i TIME, dvs. indeks beregnes i
# mellem to på hinanden værdier i TIME. Hvis TIME ikke er numerisk,
# antages at rækkefølgen er den ønskede.


malmquist <- function(X, Y, ID, TIME, 
         RTS="vrs", ORIENTATION="in",
         SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
         TRANSPOSE=FALSE, FAST=TRUE, LP=FALSE, CONTROL=NULL, LPK=NULL)  
{

   # De tidspunkter/perioder der er i data
   time <- unique(TIME)
   unit <- unique(ID)
   
   # |time| sorteres hvis variablen er numerisk ellers bruges implicit rækkefølge
   if ( is.numeric(time) )  time <- sort(time)

   # Første årstal |time| vedbliver med at være NA
   # Skal rækkefølgen ikke være som i X og Y?
   Malm <- matrix(NA, nrow=length(time), ncol=length(unit))
   TC <- matrix(NA, nrow=length(time), ncol=length(unit))
   EC <- matrix(NA, nrow=length(time), ncol=length(unit))

   E00 <- matrix(NA, nrow=length(time), ncol=length(unit))
   E01 <- matrix(NA, nrow=length(time), ncol=length(unit))
   E10 <- matrix(NA, nrow=length(time), ncol=length(unit))
   E11 <- matrix(NA, nrow=length(time), ncol=length(unit))

   # Løb perioderne igennem og beregn Malmquist for parvise perioder
   for ( t in 2:length(time) )  {
       # Find units i periode 0 og periode 1
       id0 <- ID[time[t-1]==TIME]
       id1 <- ID[time[t]==TIME]

       X0 <- X[time[t-1]==TIME,, drop=FALSE] 
       Y0 <- Y[time[t-1]==TIME,, drop=FALSE]
       X1 <- X[time[t]==TIME,, drop=FALSE]
       Y1 <- Y[time[t]==TIME,, drop=FALSE]

#print("0. periode")
#print(cbind(X0=X0, Y0=Y0, id0))
#print("1. periode")
#print(cbind(X1=X1, Y1=Y1, id1))

       m <- malmq(X0,Y0,id0, X1,Y1,id1,  RTS=RTS, ORIENTATION=ORIENTATION,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)

       # Skal rækkefølgen ikke være som i X og Y?
       Malm[t,unit %in% m$id] <- m$m
       TC[t,unit %in% m$id] <- m$tc
       EC[t,unit %in% m$id] <- m$ec

       E00[t,unit %in% m$id] <- m$e00
       E01[t,unit %in% m$id] <- m$e01
       E10[t,unit %in% m$id] <- m$e10
       E11[t,unit %in% m$id] <- m$e11
       
   }  # for (t)

   return(list(m=Malm, tc=TC, ec=EC, id=unit, time=time,
               e00=E00, e10=E10, e11=E11, e01=E01))
} # function