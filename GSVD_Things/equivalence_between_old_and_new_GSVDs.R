library(GSVD)
data(authors, package="GSVD")
Observed <- authors/sum(authors)
row.w <- rowSums(Observed)
row.W <- diag(1/row.w)
col.w <- colSums(Observed)
col.W <- diag(1/col.w)
Expected <- row.w %o% col.w
Deviations <- Observed - Expected



ca.res_gplssvd <- gplssvd(Deviations, Deviations, row.W, row.W, col.W, col.W)
# ca.res_gplssvd_2 <- gplssvd_2(Deviations, Deviations, row.W, row.W, col.W, col.W)
# ca.res_gplssvd_2 <- gplssvd_2(Deviations, Deviations, 1/row.w, 1/row.w, 1/col.w, 1/col.w)
ca.res_gplssvd_2 <- gplssvd_2(Deviations, Deviations, row.W, row.W, col.W, col.W)

ca.res_gsvd <- gsvd(Deviations,row.W,col.W)
# ca.res_gsvd_2 <- gsvd(Deviations,row.W,col.W)
# ca.res_gsvd_2 <- gsvd_2(Deviations,1/row.w,1/col.w)
ca.res_gsvd_2 <- gsvd_2(Deviations,row.W,col.W)

ca.res_geigen <- geigen(t(Deviations) %*% diag(1/row.w) %*% Deviations, col.W)
# ca.res_geigen_2 <- geigen_2(t(Deviations) %*% diag(1/row.w) %*% Deviations, col.W)
# ca.res_geigen_2 <- geigen_2(t(Deviations) %*% diag(1/row.w) %*% Deviations, 1/col.w)
ca.res_geigen_2 <- geigen_2(t(Deviations) %*% diag(1/row.w) %*% Deviations, col.W)


ca.res_geigen_alt <- geigen(Deviations %*% diag(1/col.w) %*% t(Deviations), row.W)
ca.res_geigen_alt_2 <- geigen_2(Deviations %*% diag(1/col.w) %*% t(Deviations), row.W)




mapply("/", ca.res_geigen, ca.res_geigen_2)

mapply("/", ca.res_geigen_alt, ca.res_geigen_alt_2)

mapply("/",
       ca.res_gsvd[intersect(names(ca.res_geigen_2), names(ca.res_gsvd))],
       ca.res_geigen_2[intersect(names(ca.res_geigen_2), names(ca.res_gsvd))])



mapply("/", ca.res_gplssvd, ca.res_gplssvd_2)

# mapply("/",
#        ca.res_gsvd[intersect(names(ca.res_gplssvd_2), names(ca.res_gsvd))],
#        ca.res_gplssvd_2[intersect(names(ca.res_gplssvd_2), names(ca.res_gsvd))])

mapply("/", ca.res_gsvd, ca.res_gsvd_2)

### make a new test that has a vector and then another with a real matrix, not just a diag


### make a new version of this, but test the eigen against the SVD
