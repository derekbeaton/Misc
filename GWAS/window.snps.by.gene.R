require(biomaRt)

window.snps.by.gene <- function(GET.THESE.GENES,window.size=25000){

  # snp.mart <- useEnsembl(biomart="snp", dataset="hsapiens_snp")
	# gene.mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  # snp.mart <- useMart("snp", host = "http://aug2017.archive.ensembl.org", dataset = "hsapiens_snp")
  # gene.mart <- useMart("ensembl", host = "http://aug2017.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

  snp.mart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
  gene.mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
	locations <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position', 'end_position', 'strand'), filters ='hgnc_symbol', values = GET.THESE.GENES, mart = gene.mart)
	
	ensemblids <- locations[,1]
		## I have the chr. start and end locations but I should use the latest build because rs# are static but chr. pos is not.
	rsid.mats <- vector("list",length(ensemblids))
	names(rsid.mats) <- ensemblids 
	 
	window.snps <- vector("list",length(ensemblids))	
	 
	for(i in 1:length(ensemblids)){
		
		this.gene <- getBM(attributes=c("refsnp_id","chr_name","chrom_start","chrom_end","ensembl_gene_stable_id"),filters=c("ensembl_gene"), values=c(ensemblids[i]), mart= snp.mart)
		
		
		window.snps[[i]] <- getBM(
						attributes=c('refsnp_id','chr_name','chrom_start','chrom_end'), 
						filters = c('chr_name','start','end'),
						values=list(unique(this.gene[,2]), min(this.gene[,3])-window.size, max(this.gene[,3])+window.size),
						mart=snp.mart)	
		
		#print(i)
		print(ensemblids[i])
	}
	names(window.snps) <- ensemblids
	return(window.snps)
}

