
library(ExPosition)
library(GSVD) 
	##devtools::install_github("derekbeaton/ExPosition-Family",subdir="/ExPosition2/GSVD/Package")

conditional.mca.reconstruct <- function(DATA,CONFOUND){
  
  N <- sum(DATA)
  Or <- DATA/N
  m <- rowSums(Or)
  w <- colSums(Or)
  Er <- m %o% w 
  
  MODEL <- (CONFOUND/sum(DATA)) %*% t( t(DATA) %*% CONFOUND %*% diag(1/colSums(CONFOUND)) )   
  DATA.prime <- (Or - MODEL + Er) * N
  
  ## we could just run CA here... probably will in ExPosition2
  return(DATA.prime)
}

add.alpha <- function(col, alpha=0.65){
  apply(
    sapply(col, col2rgb)/255, 
    2, 
    function(x){
      rgb(x[1], x[2], x[3], alpha=alpha)
    }
  )  
}


simple.can.ca <- function(R,row.weights,col.weights){
	
	res <- gsvd(R, row.weights%^%-1, col.weights%^%-1)
    ## first pair of vectors & value are trivial (i.e., mean), and should be removed.
    res$d.orig <- res$d.orig[-1]
    res$tau <- (res$d.orig^2 / sum(res$d.orig^2)) * 100
    res$d <- res$d[-1]
    res$fi <- res$fi[,-1]
    res$p <- res$p[,-1]
    res$u <- res$u[,-1]
    res$fii <- res$fii[,-1]
    res$fj <- res$fj[,-1]
    res$q <- res$q[,-1]
    res$v <- res$v[,-1]
    res$constraints <- minmaxHelper( rbind(res$fi, res$fj) )

	return(res)	
	
}