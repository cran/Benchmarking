# $Id: bootStat.R 72 2010-09-11 17:06:14Z Lars $

# Calculates the critical value at level |alfa| for the vector of
# trials |s|
critValue <- function(s,alfa=0.05) {
  ss_ <- sort(s)
  mean( ss_[floor(alfa*length(s))], ss_[ceiling(alfa*length(s))] )
}


# Calculate the probability of a larger value than |shat| in the vector 
# of trials |s|
typeIerror <- function(shat,s) {
  reject <- function(alfa)  {
    quantile(s,alfa,names=F) - shat
  }
  uniroot(reject,c(0,1))$root
}
