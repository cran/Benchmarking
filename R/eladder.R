# $Id: eladder.R 258 2023-08-03 12:27:03Z larso $
# Dropper succesivt den peer der giver størst stigning i eff.tallet.

eladder <- function(n, X, Y, RTS="vrs", ORIENTATION="in", 
					XREF=NULL, YREF=NULL, DIRECT=NULL, param=NULL, MAXELAD=NULL)  {

	if ( is.null(XREF) )  {
		XREF <- X
	}
	if ( is.null(YREF) )  {
		YREF <- Y
	}

	idx <- NULL
	if ( missing(MAXELAD) || is.null(MAXELAD) ) {
		MAXELAD <- dim(XREF)[1]
	} else {
		if( !is.numeric(MAXELAD) ) stop("MAXELAD must be an integer")
		MAXELAD <- min(abs(MAXELAD), dim(XREF)[1])
	}
	Xn <- X[n,,drop=FALSE]
	Yn <- Y[n,,drop=FALSE]
	# 'navne' bruges til at kunne referere til units i den oprindelige XREF
	# Efterhaanden som peers udgaar, svarer indeks i XREF ikke til de brugte
	# XREF via FRONT.IDX. Derfor bruges 'navne' til de oprindelige indeks.
	navne <- rownames(XREF)
	if ( is.null(navne) )
		navne <- rownames(YREF) 
	if ( is.null(navne) )  {
		# Hvis der ikke er nogen raekkenavne laves en loebende raekkenr.
		navne <- 1:dim(XREF)[1]
	}
	elad <- rep(NA, MAXELAD)
	for ( i in 1:MAXELAD )  {
		if ( length(idx) == MAXELAD )  {
			break  
		}
		# Brug FRONT.IDX for at kunne bruge de oprindelige indeks i X og Y
		## if (LP) cat("**", i, ": idx = ", paste(idx, collapse=", "), "\n", sep="")
		e <- dea(Xn, Yn, RTS=RTS, ORIENTATION=ORIENTATION,
				 XREF=XREF, YREF=YREF, FRONT.IDX=idx, DIRECT=DIRECT, 
				 param=param, LP=FALSE)
		# Er der nogen peers overhovedet ellers kan vi bare slutte nu
		if ( is.na(eff(e)) ) break
		if ( abs(eff(e)) == Inf ) break
		if ( is.na(peers(e)[1]) ) break
		# Gem den relevante eff.værdi
		elad[i] <- eff(e)
		# Find alle peers og find så den mest betydningsfulde
		fb <- peers(e)
		# 'fb' er indeks i XREF[idx,] og skal laves om til indeks i XREF
		## if (LP) print(fb)
		efb <- rep(NA, length(fb))
		ip <- 0
		for ( f in fb )  {
			# 'fb' og dermed 'f' er indeks i XREF[idx,] og skal laves om til indeks i XREF
			if (is.null(idx))  {
				str <- navne[f]
			} else {
				str <- navne[idx][f]
			}
			ff <- which( navne %in% str )
			e <- dea(Xn, Yn, RTS=RTS, ORIENTATION=ORIENTATION,
				 XREF=XREF, YREF=YREF, FRONT.IDX=c(idx, -ff), DIRECT=DIRECT, FAST=TRUE,
				 param=param, LP=FALSE)
			## if (LP) cat(i, ":  str = ", str, ", f = ", f, ", ff = ", ff, "; e = ", e, "\n", sep="")
			ip <- ip + 1
			efb[ip] <- e
		}
		# Peer med største effekt på eff.tal
		ip <- which.max(efb)
		# Find indeks for placering af fb[ip] i rownames(XREF)/rownames(YREF)
		if (is.null(idx))  {
			ipp <- which( navne %in% navne[fb[ip]] )
		} else {
			ipp <- which( navne %in% navne[idx][fb[ip]] )
		}		
		idx <- c(idx, -ipp)
		## if (LP) cat("*", i, ": ipp = ", ipp, ", f = ", fb[ip], "; e = ", efb[ip], "; idx = ", paste(idx, collapse=", "), "\n", sep="")
		# Hvis eff er over 1 så stopper videre beregninger; ville de give mening?
		if (efb[ip] > 1) break
	}
	elad <- elad[!is.na(elad)]
	e <- dea(Xn, Yn, RTS=RTS, ORIENTATION=ORIENTATION,
			 XREF=XREF, YREF=YREF, FRONT.IDX=idx, DIRECT=DIRECT, 
			 param=param, LP=FALSE)
	if ( is.null(idx) )  {
		idx <- NA
	} else {
		idx <- -idx
	}
	return(list(eff=elad, peer=idx, lastp=peers(e)))
}  ## eladder
