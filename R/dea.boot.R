# $Id: dea.boot.R 85 2010-11-06 22:18:27Z Lars $

# Bootstrap DEA functions, a wrapper for FEAR::boot.sw98


# boot.sw98(XOBS, YOBS, NREP = 2000, DHAT = NULL, 
#    RTS = 1, ORIENTATION = 1, alpha = 0.05, CI.TYPE=2,
# 	  XREF = NULL, YREF = NULL, DREF = NULL, 
# 	  OUTPUT.FARRELL = FALSE, NOPRINT = FALSE, errchk = TRUE)





dea.boot <- function(X,Y, NREP=50, EFF=NULL, RTS="vrs", ORIENTATION="in",
       alpha=0.05, XREF=NULL, YREF=NULL, EREF=NULL)  
{

   if ( !("FEAR" %in% .packages(T)) )
      stop("Note: dea.boot only works if the package FEAR is installed")


   rts_ <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( is.real(RTS) )  {
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
   }
   RTS <- tolower(RTS)
   rts <- which(RTS == rts_) -1
   if ( rts < 1 || 3 < rts )
      stop("Invalid value of RTS in call to dea.boot")


   orientation_ <- c("in-out","in","out","graph")
   if ( is.real(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   orientation <- which(ORIENTATION == orientation_) -1
   if ( orientation < 1 || 3 < orientation )
      stop("Invalid value of ORIENTATION in call to dea.boot")

   if (ORIENTATION=="out") farrell<-TRUE else farrell<-FALSE

   # print(paste("rts =",rts,"   orientation=",orientation))

   b <- NA
   tryCatch( b <- FEAR::boot.sw98(t(X), t(Y), NREP, EFF, rts, 
                    orientation, alpha,,XREF, YREF, EREF,
                    OUTPUT.FARREL=farrell) , 
       warning = function(w) print(w),
       error = function(e) {
              print(e)
              stop("dea.boot aborted:  FEAR not installed or loaded")
       },
       finaly = print("FEAR::boot.sw98 finished with bootstrap",quote=FALSE)
   )

   # boot.sw98(t(x), t(y), NREP=NREP)

   if ( farrell )  {
      # bb <- b
      bb <- list(bias=b$bias, var=b$var, conf.int=b$conf.int, 
               eff=b$dhat, eff.bc=b$dhat.bc, eref=b$dref, boot=b$boot)
  } else {
      # print("Omregner til Farrell")
      # Omregn til Farrell efficiencer når det nu ikke er Farrell
      eff <- 1/b$dhat
      boot <- 1/b$boot 
      bias <- rowMeans(boot) - eff
      eff.bc <- eff - bias

      boot_ <- sweep(boot,1,FUN="-",eff)
      ci_ <- 
        t(apply(boot_,1, quantile, 
          probs=c(1-0.5*alpha, 0.5*alpha), type=9, na.rm=TRUE))
      ci <- sweep( -ci_, 1, FUN="+", eff )
      rm(boot_, ci_)
 
      var <- apply(boot,1,var)
      if ( is.null(b$dref) ) eref <- NULL  else  eref <- 1/b$dref
      bb <- list(bias=bias, var=var, conf.int=ci, 
                 eff=eff, eff.bc=eff.bc, eref=eref, boot=boot)
   }

   return( bb )
}

