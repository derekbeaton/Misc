## pca from GPLSSVD

library(GSVD)
data("snps.druguse", package="GSVD")

model.matrix(~ ., data=cbind(snps.druguse$DATA1,snps.druguse$DATA2),
             contrasts.arg = lapply(cbind(snps.druguse$DATA1,snps.druguse$DATA2), contrasts, contrasts=FALSE))[,-1]


 data("snps.druguse", package="GSVD")
 X <- model.matrix(~ ., data=snps.druguse$DATA1, contrasts.arg = lapply(snps.druguse$DATA1, contrasts, contrasts=FALSE))[,-1]
 Observed <- X/sum(X)
 row.w <- rowSums(Observed)
 col.w <- colSums(Observed)
 Expected <- row.w %o% col.w
 Deviations <- Observed - Expected
 ca.res <- gsvd(Deviations,1/row.w,1/col.w)


epres <-  ExPosition::epMCA(snps.druguse$DATA1, graphs=F, correction = NULL)
epres1 <-  ExPosition::epCA(X, graphs=F)



epres <-  ExPosition::epMCA( cbind(snps.druguse$DATA1,snps.druguse$DATA2), graphs=F)
epres1 <- ExPosition::epCA(ExPosition::makeNominalData(cbind(snps.druguse$DATA1,snps.druguse$DATA2)), graphs=F)

