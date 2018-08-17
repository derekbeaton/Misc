## these go along with the appendix (part 1)

load('tcata.sim.data.rda')
source('helpers.R')

### This file helps clarify the first trivial (mean) vector from CanCA by comparison with standard CA.

Y <- tcata.sim.data$DATA
X.products <- makeNominalData(as.matrix(tcata.sim.data$DESIGNS[,1]))
X.time <- makeNominalData(as.matrix(tcata.sim.data$DESIGNS[,2]))
X.products.time <- makeNominalData(tcata.sim.data$DESIGNS)



### Appendix A.1: Centered and uncentered CA and CanCA
## standard CA (through the GSVD)
Oy <- Y / sum(Y)
wi <- rowSums(Oy)
	WI <- diag(wi)
wj <- colSums(Oy)
	WJ <- diag(wj)
Ey <- wi %o% wj
Zy <- Oy - Ey
standard.ca <- gsvd(Zy, WI %^% -1, WJ %^% -1)



## uncentered standard CA, where the first principal axis is trivial.
WI.Y <- diag(rowSums(Y))
WJ.Y <- diag(colSums(Y))
standard.ca_uncentered <- gsvd(Oy, WI %^% -1, WJ %^% -1)

	## a quick look at the singular values
standard.ca$d.orig
standard.ca_uncentered$d.orig
	## a quick look at column scores
standard.ca$fj
standard.ca_uncentered$fj


	## reconstruct the expected values from the uncentered results:
Ey_from.uncentered <- (standard.ca_uncentered$p[,1] %*% t(standard.ca_uncentered$q[,1])) 
Ey / ( Ey_from.uncentered)
	## this shows that the first pair of vectors creates the expected matrix.

standard.ca_mixture <- gsvd(Oy - Ey_from.uncentered, WI %^% -1, WJ %^% -1)


## Canonical CA -- the "uncentered" version (as presented in the paper)
	## I need to get this to be exactly the same as the _alt below; then re-do the equations a bit. I should be close
WI.star <- t(X.products.time) %*% WI %*% X.products.time
	WI.star <- WI.star / sum(WI.star)
R <- t(X.products.time) %*% Y
	Or <- R / sum(R)
can.ca.res <- gsvd(Or, WI.star %^% -1, WJ %^% -1)

## We can now build the expected matrix for CanCA
Er <- (can.ca.res$p[,1] %*% t(can.ca.res$q[,1]))	
## and thus the deviations for CanCA
Zr <- Or - Er
## and then apply that through the GSVD for a "centered" version:
can.ca.res_centered <- gsvd(Zr, WI.star %^% -1, WJ %^% -1)


	## a quick look at the singular values
can.ca.res$d.orig
can.ca.res_centered$d.orig
	## a quick look at column scores
can.ca.res$fj
can.ca.res_centered$fj





### Appendix A.2: Component-wise contributions to the variance
	### only illustrated here with CA

## From the GSVD code here, P & Q are the generalized vectors, whereas U & V are the standard vectors from the svd().

deweighted.gen.vector <- ((WJ %^% (-1/2)) %*% standard.ca$q)
	# show equivalence to standard vectors:
deweighted.gen.vector / standard.ca$v

contributions <- deweighted.gen.vector * deweighted.gen.vector