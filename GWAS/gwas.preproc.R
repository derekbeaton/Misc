
gwas.preproc <- function(data.in,geno=0.1,mind=0.1,maf=0.05){
	
	##zero-th: find any columns that are only 1 genotype.
	na.ret.thing <- apply(data.in,2,function(x){unique(x[!is.na(x)])})
		if(is.list(na.ret.thing)){
		single.genos <- which( unlist( lapply(na.ret.thing,length) )<2 )
		if(length(single.genos)>=1){	
			data.in <- data.in[,-c(single.genos)]
		}
	}
	##first: SNP missingness
	geno.missing <- which( colSums(is.na(data.in)) > nrow(data.in)*geno )
	if(length(geno.missing)>=1){	
		data.in <- data.in[,-c(geno.missing)]
	}

	##second: person missingness
	mind.missing <- which( rowSums(is.na(data.in)) > ncol(data.in)*mind )
	if(length(mind.missing)>=1){
		data.in <- data.in[-c(mind.missing),]
	}
		
	##third: minor allele frequencies
	ma <- tapply(data.in,col(data.in),function(x){table(unlist(strsplit(as.vector(x),"")))},simplify=F)	
	names(ma) <- colnames(data.in)	
	mafreqs <- lapply(ma,function(x){x/sum(x)})
	maf.drop <- which(lapply(mafreqs,function(x){sum(x<maf)})==1)
	if(length(maf.drop)>=1){
		data.in <- data.in[,-c(maf.drop)]
	}

	  ## remove this as its not standard. will make custom code to handle this.
	# ##fourth: find new minor alleles, so we can combine genotypes.
	# ma <- tapply(data.in,col(data.in),function(x){table(unlist(strsplit(as.vector(x),"")))},simplify=F)	
	# names(ma) <- colnames(data.in)		
	# ma.finder <- lapply(ma,function(x){ a<-x/sum(x); a==min(a) }) ## apparently not used...
	# geno.counts <- tapply(data.in,col(data.in),function(x){ summary(as.factor(x[!is.na(x)])) },simplify=F)
	# 	names(geno.counts) <- colnames(data.in)
	# find.low.genos <- tapply(data.in,col(data.in),function(x){ summary(as.factor(x[!is.na(x)])) < (length(x)*maf) },simplify=F)
	# 	names(find.low.genos) <- colnames(data.in)
	# 
	# find.combos <- lapply(find.low.genos,function(x){
	# 									a=which(x); ##tells me which item is below threshold
	# 									if(length(a)==0){
	# 										thing=NULL
	# 									}else{
	# 										thing=unique(names(x)[c(2,a)]) ##tells me the two items to combine into one...?
	# 									}
	# 								})
	# 								
	# change.these <- find.combos[which(unlist(lapply(find.combos,function(x){!is.null(x)})))]
	# for(i in 1:length(change.these)){
	# 	this.one <- change.these[i]
	# 	change.to <- paste(unlist(this.one),collapse="+")
	# 	data.in[c(which(data.in[,names(this.one)]==unlist(this.one)[1]),which(data.in[,names(this.one)]== unlist(this.one)[2])),names(this.one)] <- change.to		
	# }
	
	return(data.in)
}


