
Observed <- authors/sum(authors)
row.w <- rowSums(Observed)
row.W <- diag(1/row.w)
col.w <- colSums(Observed)
col.W <- diag(1/col.w)
Expected <- row.w %o% col.w
Deviations <- Observed - Expected
ca.res_gsvd <- gsvd(Deviations,row.W,col.W)
ca.res_geigen <- geigen(t(Deviations) %*% diag(1/row.w) %*% Deviations, col.W)

### quick multiplication benchmraks

test_mat <- t(Deviations) %*% diag(1/row.w) %*% Deviations
rownames(test_mat) <- paste0("ROW: ", rownames(test_mat))
colnames(test_mat) <- paste0("COL: ", colnames(test_mat))

t(Deviations / row.w) %*% Deviations



DAT <- sweep(sweep(test_mat, 2, sqrt(1/col.w), "*"), 1, sqrt(1/col.w), "*")



sweep(sweep(test_mat, 2, sqrt(col.w), "*"), 1, sqrt(col.w), "*")


N <- 1000
P <- 100
DD <- matrix(rnorm(N * P), N, P)
RW <- rnorm(N)
mbm_out <- microbenchmark::microbenchmark(

  element_mult = t(DD / RW),
  matrix_mult = t(DD) %*% diag(1/RW),
  times = 1000

)

 t(test_mat / sqrt(col.w)) / sqrt(col.w)


DDD <- matrix(rnorm(P * P), P, P)
CW <- rnorm(P)
CW <- CW + (abs(min(CW))+1)

A <- sweep(sweep(DDD, 2, sqrt(1/CW), "*"), 1, sqrt(1/CW), "*")
B <- t(t(DDD / sqrt(CW)) / sqrt(CW))
C <- diag(1/sqrt(CW)) %*% DDD %*% diag(1/sqrt(CW))

mbm_out <- microbenchmark::microbenchmark(

  sweep_the_leg = sweep(sweep(DDD, 2, sqrt(1/CW), "*"), 1, sqrt(1/CW), "*"),
  vect_1 = t(t(DDD / sqrt(CW)) / sqrt(CW)),
  matrix_mult = diag(1/sqrt(CW)) %*% DDD %*% diag(1/sqrt(CW)),
  times = 1000

)



N <- 100
P <- 1000
DD <- matrix(rnorm(N * P), N, P)
mbm_out <- microbenchmark::microbenchmark(
  swept_up = {svd_res <- svd(DD); sweep(svd_res$u, 2, svd_res$d, "*")},
  vectorz = {svd_res <- svd(DD); t((t(svd_res$u) * svd_res$d))},
  matriz = {svd_res <- svd(DD); svd_res$u %*% diag(svd_res$d)},
  times = 1000

)
