#PURPOSE: Get lists of genes associated with each germ cell cluster in Green et al. 2018 (mouse spermatogenesis single cell RNAseq paper); convert to protein lists

library(biomaRt)

#Read in data
mydata<-read.csv("hermannEtal2018_1-s2.0-S2211124718316024-mmc2.csv", header=TRUE, comment.char="#")
print(head(mydata))

#Load biomaRt
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('external_gene_name','ensembl_peptide_id'), mart=ens_mus)
colnames(all_genes)<-c("protein_id", "gene_name")
print(head(all_genes))

myGCs<-unique(mydata$cluster)
print(myGCs)
#Also assign clusters to cell types based on Hermann et al. Fig. 1E
spg<-c()
prelep<-c()
lz_pac<-c()
pac_dip<-c()
dip_rs<-c()
rs<-c()
sert<-c()
peri<-c()
for(i in myGCs){
	geneNames<-mydata$gene[which(mydata$cluster==i)]
	protIDs<-all_genes$protein_id[which(all_genes$gene_name %in% geneNames)]
	protIDs_unique<-unique(protIDs)
	protIDs_filtered<-protIDs_unique[which(protIDs_unique != "")]
	print(paste("Number of unique protein IDs for germ cell cluster", i, length(protIDs_filtered)))
	write(protIDs_filtered, file=paste0("prot_list.GCcluster",i,".txt"), ncolumns=1)
	if(i == 16){
		spg<-c(spg, protIDs_filtered)
	} else if(i == 10){
		prelep<-c(prelep, protIDs_filtered)
	} else if(i == 8){
		lz_pac<-c(lz_pac, protIDs_filtered)
	} else if(i == 6){
		pac_dip<-c(pac_dip, protIDs_filtered)
	} else if(i == 9){
		dip_rs<-c(dip_rs, protIDs_filtered)
	} else if(i %in% c(12, 2, 7, 13, 11, 5, 3, 1, 4)){
		rs<-c(rs, protIDs_filtered)
	} else if(i %in% c(14, 15)){
		sert<-c(sert, protIDs_filtered)
	} else if(i == 17){
		peri<-c(peri, protIDs_filtered)
	} else{
		print(paste("Warning: no cell type for cluster", i))
	}
}

write(unique(spg), file="prot_list.hermannEtal2018.spermatogonia.txt", ncolumns=1)
write(unique(prelep), file="prot_list.hermannEtal2018.prelep.txt", ncolumns=1)
write(unique(lz_pac), file="prot_list.hermannEtal2018.lz_pac.txt", ncolumns=1)
write(unique(pac_dip), file="prot_list.hermannEtal2018.pac_dip.txt", ncolumns=1)
write(unique(dip_rs), file="prot_list.hermannEtal2018.dip_rs.txt", ncolumns=1)
write(unique(rs), file="prot_list.hermannEtal2018.rs.txt", ncolumns=1)
write(unique(sert), file="prot_list.hermannEtal2018.sert.txt", ncolumns=1)
write(unique(peri), file="prot_list.hermannEtal2018.peri.txt", ncolumns=1)

print("Done with 05_getHermannEtal_GCcluster_proteins.r")
