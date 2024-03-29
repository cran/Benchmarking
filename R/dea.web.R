# $Id: dea.web.R 244 2022-05-05 14:31:31Z X052717 $

# Foerste forsoeg med en web-graf for efficiens og tilhoerende input/output

dea.web <-
function(X, E, N=NULL, txt=NULL, add=FALSE, 
            GRID=FALSE, fex=1,
            # RANGE=FALSE, param=NULL, 
            ..., xlim)
# x er input eller output for en firm eller for firm n
#
{

   if (dim(X)[1] > 1 && missing(N) )
      stop("If X is for more firms then N must be present")
   if (dim(X)[1] > 1 && missing(N) )
      stop("If X is for more firms then N must be present")

   sl <- 0
   if ( is(E, "Farrell") )  {
      # Hvordan skelne mellem sx eller st der er behov for?
    if ( !is.null(E$sx) && E$ORIENTATION=="in" )  sl <- E$sx
    else if ( !is.null(E$sx) && E$ORIENTATION=="out" )  sl <- -E$sy
       E <- eff(E)
   }
   if ( length(E) != dim(X)[1] )
      stop("Firms in X and E must be the same")

   if ( is(X, "matrix") )  {
      m <- dim(X)[2]
      X <- X[N,]
      E <- E[N]
   } else  m <- length(X)
   # Nu er en X et array
   if ( is(sl, "matrix") )  {
      sl <- sl[N,]
   }

   step <- 2*pi/m
   angle <- seq(from=0, to=2*pi-step, by=step)

   dots = list(...)
   if ( missing(xlim) )  xlim <- ylim <- 
                      max(1,E) * (c(-max(X), max(X))+.01) * 1.1

   plot( c(X * cos(angle), X[1]), c(X * sin(angle),0), type="l",
        xlim=xlim, ylim=ylim)
   lines( c((E * X-sl) * cos(angle), E*X[1]-sl[1]), c((E * X-sl) * sin(angle),0))

   # Polaere linjer for hver vare
   segments(rep(0,m),rep(0,m), xlim[2]*cos(angle), xlim[2]*sin(angle), 
    col="darkgray")

   text(xlim[2]*cos(angle)*1.05, xlim[2]*sin(angle)*1.05, 1:m) #, col="darkgray")

  if ( GRID )  {
       grid(col="darkgray")
       box(col="grey")
  }
  if ( is(txt, "logical") && txt )  {
      if ( !is(X, "matrix") )  {
         if ( !is.null(rownames(X)) )  {
            txt <- rownames(X)
         } else {
            txt <- 1:dim(X)[1]
         }
      } else {
         txt <- 1:length(X)
      }
   }
   if ( !is(txt, "logical") && length(txt) > 0 ) {
     # Evt. tekst paa punkter saettes lidt nede til hoejre
     text(X,X,txt,adj=c(-.75,.75),cex=fex)
   }

}  # dea.web

