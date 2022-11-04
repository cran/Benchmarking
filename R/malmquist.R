# $Id: malmquist.R 247 2022-09-15 22:38:27Z X052717 $

# Beregner Malmquist indeks for enhederne ID over tidspunkterne i TIME

# Det forudsaettes at der ikke er huller i TIME, dvs. indeks beregnes i
# mellem to paa hinanden vaerdier i TIME. Hvis TIME ikke er numerisk,
# antages at raekkefoelgen er den oenskede.

# Parvise perioder beregnes via functionen 'malmq'


malmquist <- function(X, Y, ID, TIME, 
         RTS="vrs", ORIENTATION="in", SAMEREF=FALSE,
         SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
         TRANSPOSE=FALSE, FAST=TRUE, LP=FALSE, CONTROL=NULL, LPK=NULL)  
{

   # De tidspunkter/perioder der er i data
   time <- unique(TIME)
   unit <- unique(ID)
   
   # |time| sorteres hvis variablen er numerisk ellers bruges implicit raekkefoelge
   # Boer saa ikke al data sorteres? Nej fordi det er ikke TIME der sorteres, men alene
   # de entydige TIME vaerdier uafhaenhaengigt at units.
   if ( is.numeric(time) )  time <- sort(time)

   # Vaerdier for foerste aarstal |time| vedbliver med at vaere NA
   # Raekkefoelgen er som i ID?
   Malm <- array(NA, dim=c(length(ID)))
   TC <- array(NA, dim=c(length(ID)))
   EC <- array(NA, dim=c(length(ID)))

   E00 <- array(NA, dim=c(length(ID)))
   E01 <- array(NA, dim=c(length(ID)))
   E10 <- array(NA, dim=c(length(ID)))
   E11 <- array(NA, dim=c(length(ID)))

   # Loeb perioderne igennem og beregn Malmquist for parvise perioder
   for ( t in 2:length(time) )  {
        cat("Period ",t,"\n")
        flush.console()
       # Find units i periode 0 og periode 1
       id0 <- ID[time[t-1]==TIME]
       id1 <- ID[time[t]==TIME]

       X0 <- X[time[t-1]==TIME,, drop=FALSE] 
       Y0 <- Y[time[t-1]==TIME,, drop=FALSE]
       X1 <- X[time[t]==TIME,, drop=FALSE]
       Y1 <- Y[time[t]==TIME,, drop=FALSE]
       
       m <- malmq(X0,Y0,id0, X1,Y1,id1,  RTS=RTS, ORIENTATION=ORIENTATION,
         SAMEREF=SAMEREF, SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)

    # Raekkefoelgen skal vaere som i X og Y, dvs. som ID
    for (i in m$id)  {
        Malm[ID==i & TIME==time[t]] <- m$m[m$id==i] # Malmquist indeks for aendring i produktivitet
        TC[ID==i & TIME==time[t]] <- m$tc[m$id==i] # teknisk aendring, flytning af frontier
        EC[ID==i & TIME==time[t]] <- m$ec[m$id==i] # aendring i effektivitet

        E00[ID==i & TIME==time[t]] <- m$e00[m$id==i]
        E01[ID==i & TIME==time[t]] <- m$e01[m$id==i]
        E10[ID==i & TIME==time[t]] <- m$e10[m$id==i]
        E11[ID==i & TIME==time[t]] <- m$e11[m$id==i]
        }

        if ( t==2)  {
            # Foerste aar kan saettes 
            for (i in m$id)  {
                E11[ID==i & TIME==time[t-1]] <- m$e00[m$id==i]
            }
        }
    }  # for (t)
   

   return(list(m=Malm, tc=TC, ec=EC, id=ID, time=TIME,
               e00=E00, e10=E10, e11=E11, e01=E01))
} # function
