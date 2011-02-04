# $Id: bootStat.R 102 2011-01-23 18:50:21Z Lars $

# Calculates the critical value at level |alpha| for the vector of
# trials |s|
critValue <- function(s, alpha=0.05) {
  if ( alpha <= 0 || alpha >= 1 ) 
     stop("The argument alpha must be between 0 and 1")
  ss_ <- sort(s)
  mean( ss_[floor(alpha*length(s))], ss_[ceiling(alpha*length(s))] )
}


# Calculate the probability of a larger value than |shat| in the vector 
# of trials |s|
typeIerror <- function(shat,s) {
  reject <- function(alfa)  {
    quantile(s, alfa, names=F) - shat
  }
  if ( reject(0) * reject(1) > 0 )  {
      # Ingen loesning til unitroot, saa enten 0% eller 100%
      if ( shat <= min(s) )  return(0)
      if ( shat >= max(s) )  return(1)
  }
  uniroot(reject,c(0,1))$root
}
