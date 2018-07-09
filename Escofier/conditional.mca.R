	## via Brigitte Escofier
conditional.mca <- function(X,Y){
	####FOR GLOBAL USE
	#X <- makeNominalData(dat)
	N <- sum(X)
	m <- rowSums(X)/N
	w <- colSums(X)/N
	Or <- X/N
	Er <- m %o% w
	DATA <- Or
	MARGINS <- Er
	####FOR GLOBAL USE
	
	C <- t(X) %*% Y
	kt <- colSums(Y)
	ktmat <- matrix(kt,nrow(C),ncol(C),byrow=T)
	pre.model <- ((Y/N) %*% t(C/ktmat))

	MODEL <- pre.model
	X.prime <- (DATA - MODEL + MARGINS) * N
	return(X.prime)
}
