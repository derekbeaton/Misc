load('tcata.sim.data.rda')
source('helpers.R')

# Y <- DATA
# X.time <- makeNominalData(TIME)
# X.products <- makeNominalData(PRODUCTS)
# X.products.time <- makeNominalData(PRODUCTS.and.TIME)


Y <- tcata.sim.data$DATA
X.products <- makeNominalData(as.matrix(tcata.sim.data$DESIGNS[,1]))
X.time <- makeNominalData(as.matrix(tcata.sim.data$DESIGNS[,2]))
X.products.time <- makeNominalData(tcata.sim.data$DESIGNS)


attr.cols <- 1#"grey40"
prod.cols <- c("#023FA5", "#BB7784", "#E28912", "#0FCFC0")
products <- c("PA", "PB", "PC", "PD")



## standard CA
standard.ca.res <- epCA(Y,DESIGN=X.products,make_design_nominal = F,graphs=F)

## Canonical CA
WI <- t(X.products.time) %*% diag(rowSums(Y)) %*% X.products.time
WJ <- diag(colSums(Y))
can.ca.res <- simple.can.ca(t(X.products.time) %*% Y,WI,WJ)


## Conditional CA
    ## step 1: remove effect via Escofier method
Y.prime <- conditional.mca.reconstruct(Y,X.time)
    ## step 2: CA
con.ca.res <- epCA(Y.prime,DESIGN=X.products,make_design_nominal = F,graphs=F)  



standard.ca.dots <- 1.5 * (rowSums(standard.ca.res$ExPosition.Data$cj[,1:2]) - min(rowSums(standard.ca.res$ExPosition.Data$cj[,1:2]))) / (max(rowSums(standard.ca.res$ExPosition.Data$cj[,1:2])) - min(rowSums(standard.ca.res$ExPosition.Data$cj[,1:2]))) + 1

can.ca.dots <- 1.5 * (rowSums((can.ca.res$v^2)[,1:2]) - min(rowSums((can.ca.res$v^2)[,1:2]))) / (max(rowSums((can.ca.res$v^2)[,1:2])) - min(rowSums((can.ca.res$v^2)[,1:2]))) + 1

con.ca.dots <- 1.5 * (rowSums(con.ca.res$ExPosition.Data$cj[,1:2]) - min(rowSums(con.ca.res$ExPosition.Data$cj[,1:2]))) / (max(rowSums(con.ca.res$ExPosition.Data$cj[,1:2])) - min(rowSums(con.ca.res$ExPosition.Data$cj[,1:2]))) + 1
