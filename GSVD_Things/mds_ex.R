THIS_MAT <- as.matrix(dist(wine$objective)) #as.matrix(eurodist)


# rdc <- .Call(stats:::C_DoubleCentre, as.matrix(jocn.2005.fmri$images$data)^2)

cmd_res <- cmdscale(THIS_MAT, eig = T)

D <- as.matrix(THIS_MAT)^2
masses <- rep(1/nrow(D), nrow(D))
Mrepmat <- matrix(masses,nrow=nrow(D),ncol=ncol(D))
DATA_dimensions <- dim(D)

if(is.null(dim(masses))){
  BigXi <- diag(DATA_dimensions[1]) - (matrix(1,DATA_dimensions[1],1) %*% masses)
}else{
  BigXi <- diag(DATA_dimensions[1]) - (matrix(1,DATA_dimensions[1],1) %*% diag(masses))
}
S <- -.5 * sqrt(Mrepmat) * BigXi %*% D %*% t(BigXi) * sqrt(t(Mrepmat))

s_eigen <- eigen(S)

mds_res <- ExPosition::epMDS(D, T, graphs = F)



D <- as.matrix(dist(wine$objective))^2
masses <- rep(1/nrow(D), nrow(D))
Xi <- matrix(-masses, length(masses), length(masses))
diag(Xi) <- (1-masses)



BigXi %*% (-D) %*% t(BigXi)
geigen_res <- GSVD::geigen((-D / (nrow(D) * 2)), BigXi)

