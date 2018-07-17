source('analyses.R')

begin.locs <- seq(1,nrow(standard.ca.res$ExPosition.Data$fi),nrow(standard.ca.res$ExPosition.Data$fi)/4)
end.locs <- c(begin.locs[2:length(begin.locs)]-1,nrow(standard.ca.res$ExPosition.Data$fi))


## standard CA
prettyPlot(standard.ca.res$ExPosition.Data$fj,constraints=lapply(standard.ca.res$Plotting.Data$constraints,"*",1.2),col=attr.cols,dev.new=T,main=NULL,xlab=paste0("Component 1 (",format(round(standard.ca.res$ExPosition.Data$t[1],digits=1),nsmall=1),"%)"), ylab=paste0("Component 2 (",format(round(standard.ca.res$ExPosition.Data$t[2],digits=1),nsmall=1),"%)"), fg.col="grey", bg.col="grey", display_names=FALSE,cex=standard.ca.dots)

color <- 0
for (product in products){
    color <- color + 1
    xy <- standard.ca.res$ExPosition.Data$fi[product == tcata.sim.data$DESIGNS[,1],1:2]

    for(i in 2:nrow(xy)) {
        lines(xy[c(i-1,i),], col=prod.cols[color], lwd=3, pch=19)
    }
}

text(standard.ca.res$ExPosition.Data$fi[c(begin.locs[1],end.locs[1]),],c("PA.T1","PA.T19"),col=prod.cols[1],pos=2, font=2)
text(standard.ca.res$ExPosition.Data$fi[c(begin.locs[2],end.locs[2]),],c("PB.T1","PB.T19"),col=prod.cols[2],pos=1, font=2)
text(standard.ca.res$ExPosition.Data$fi[c(begin.locs[3],end.locs[3]),],c("PC.T1","PC.T19"),col=prod.cols[3],pos=4, font=2)
text(standard.ca.res$ExPosition.Data$fi[c(begin.locs[4],end.locs[4]),],c("PD.T1","PD.T19"),col=prod.cols[4],adj=c(-0,-0.5), font=2)

text(standard.ca.res$ExPosition.Data$fj[-4,], labels=rownames(standard.ca.res$ExPosition.Data$fj)[-4], pos=c(3,1,3,3,3), font=2)
text(t(standard.ca.res$ExPosition.Data$fj[4,1:2]), labels=rownames(standard.ca.res$ExPosition.Data$fj)[4], adj=c(0,1.5), font=2)


## canonical CA
prettyPlot(can.ca.res$fj,dev.new=T,col=attr.cols,main=NULL,constraints = lapply(can.ca.res$constraints,"*",1.2),xlab=paste0("Component 1 (",round(can.ca.res$tau,digits=1)[1],"%)"), ylab=paste0("Component 2 (",round(can.ca.res$tau,digits=1)[2],"%)"),fg.col="grey", bg.col="grey", display_names=FALSE,cex= can.ca.dots)

lines(can.ca.res$fi[grep("^TIME\\.",rownames(can.ca.res$fi)),], col="blue", lwd=3, pch=19)
text(can.ca.res$fi[c("TIME.POINT.T1","TIME.POINT.T19"),],c("T.1","T.19"),col="blue",pos=1, font=2, cex=1.5)
points(can.ca.res$fi[grep("^PRODUCT",rownames(can.ca.res$fi)),],col=prod.cols,cex=2,pch=20)
text(can.ca.res$fi[grep("^PRODUCT",rownames(can.ca.res$fi)),],c("P.A","P.B","P.C","P.D"),col= prod.cols,pos=c(1,3,3,3), font=2, cex=1.5)

text(can.ca.res$fj, labels=rownames(can.ca.res$fj), pos=c(3,1,3,1,3,3), font=2)


## conditional CA
prettyPlot(con.ca.res$ExPosition.Data$fj,col=attr.cols,dev.new=T,constraints = lapply(con.ca.res$Plotting.Data$constraints,"*",1.1),main=NULL, xlab=paste0("Component 1 (",round(con.ca.res$ExPosition.Data$t,digits=1)[1],"%)"), ylab=paste0("Component 2 (",round(con.ca.res$ExPosition.Data$t,digits=1)[2],"%)"), fg.col="grey", bg.col="grey", display_names=FALSE,cex= con.ca.dots)


color <- 0
for (product in products){
    color <- color + 1
    xy <- con.ca.res$ExPosition.Data$fi[product == tcata.sim.data$DESIGNS[,1],1:2]

    for(i in 2:nrow(xy)) {
        lines(xy[c(i-1,i),], col=prod.cols[color], lwd=3, pch=19)
    }
}

          
text(con.ca.res$ExPosition.Data$fi[c(begin.locs[1],end.locs[1]),],c("PA.T1","PA.T19"),col=prod.cols[1],pos=2, font=2)
text(con.ca.res$ExPosition.Data$fi[c(begin.locs[2],end.locs[2]),],c("PB.T1","PB.T19"),col=prod.cols[2],pos=2, font=2)
text(con.ca.res$ExPosition.Data$fi[c(begin.locs[3],end.locs[3]),],c("PC.T1","PC.T19"),col=prod.cols[3],pos=c(2,3), font=2)
text(con.ca.res$ExPosition.Data$fi[c(begin.locs[4],end.locs[4]),],c("PD.T1","PD.T19"),col=prod.cols[4],pos=c(2,4), font=2)

text(con.ca.res$ExPosition.Data$fj, labels=rownames(con.ca.res$ExPosition.Data$fj), pos=c(1,4,1,4,1,4), font=2)
