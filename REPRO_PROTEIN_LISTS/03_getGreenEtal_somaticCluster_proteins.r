#PURPOSE: Get lists of genes associated with each germ cell cluster in Green et al. 2018 (mouse spermatogenesis single cell RNAseq paper); convert to protein lists

library(biomaRt)

#Read in data
mydata<-read.csv("greenEtal2018_1-s2.0-S1534580718306361-mmc5.csv", header=TRUE, comment.char="#")
print(head(mydata))

#Load biomaRt
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('external_gene_name','ensembl_peptide_id'), mart=ens_mus)
colnames(all_genes)<-c("protein_id", "gene_name")
print(head(all_genes))

mySCs<-unique(mydata$SomaticCellType)
print(mySCs)
for(i in mySCs){
	geneNames<-mydata$gene[which(mydata$SomaticCellType==i)]
	protIDs<-all_genes$protein_id[which(all_genes$gene_name %in% geneNames)]
	protIDs_unique<-unique(protIDs)
	protIDs_filtered<-protIDs_unique[which(protIDs_unique != "")]
	print(paste("Number of unique protein IDs for somatic cell type", i, length(protIDs_filtered)))
	write(protIDs_filtered, file=paste0("prot_list.Somatic",i,".txt"), ncolumns=1)
}

print("Done with 03_getGreenEtal_somaticCluster_proteins.r")
