	##reconstruction function.
	## requires ExPosition input for now.
	
ca.recon <- function(DATA,res,comps=NULL){
	
	if(is.null(comps)){
		comps <- 1:length(res$ExPosition.Data$pdq$Dv)
	}
	if(length(comps)==1){
		comps <- 1:comps
	}
	if(length(comps)==1){
		comps <- 1:length(res$ExPosition.Data$pdq$Dv)		
	}
	
		## all of this also provides a sense of "generalized CA" as proposed by Escofier. Just an interesting note for later...
	mat1 <- res$ExPosition.Data$pdq$p[,comps] %*% diag(res$ExPosition.Data$pdq$Dv[comps]) %*% t(res$ExPosition.Data$pdq$q[,comps])

	mat2 <- matrix(res$ExPosition.Data$c,nrow(mat1),ncol(mat1),byrow=T)
	
	mat3 <- diag(res$ExPosition.Data$M * sum(DATA)) %*% (mat1+mat2)
	
	return(mat3)
}
