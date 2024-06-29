#PURPOSE: Extract lists of tissue-enriched genes from Dean et al. 2009 (MBE)

library(biomaRt)

#Read in data
f1<-read.csv("LiEtal2017_mouse_tissue_specific.csv", header=TRUE, comment.char="#")

print(paste("Number of tissue-specific genes:", nrow(f1)))
print("By tissue:")
for(i in unique(f1$tissue)){
	print(i)
	print(length(which(f1$tissue==i)))
}

#Load biomaRt
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','ensembl_peptide_id'), mart=ens_mus)
colnames(all_genes)<-c("gene_id", "protein_id")
print(head(all_genes))

#Loop through all genes and get corresponding protein ID
PIs<-c()
for(i in 1:nrow(f1)){
	#Get gene lists by tissue
	geneID<-f1$gene_id[i]
	print(geneID)
	#Change gene ID to protein IDs
	protIDs<-all_genes$protein_id[which(all_genes$gene_id == geneID)]
	#print(protIDs)
	#Get unique protein IDs only
	protIDs_unique<-unique(protIDs)
	#print(protIDs_unique)
	#Remove ""
	protIDs_filter<-protIDs_unique[which(protIDs_unique != "")]
	#If no protID, add NA
	if(length(protIDs_filter)==0){
		#print("NULL")
		PIs<-c(PIs, NA)
	} else if(length(protIDs_filter)>1){
		#print("MULT")
		protIDs_final<-paste(protIDs_filter, collapse=";")
		#print(protIDs_final)
		PIs<-c(PIs, protIDs_final)
	} else{
		#print("ONE")
		PIs<-c(PIs, protIDs_filter)
	}
	print(PIs)
}

#Append
finaldf<-as.data.frame(cbind(f1, protID=PIs))

#Write output
write.csv(finaldf, file="LiEtal2017_mouse_tissue_specific.withProtID.csv", row.names=FALSE, quote=FALSE)
