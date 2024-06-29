#PURPOSE: Test for GO enrichment for genes under selection associated with small relative testes mass

library(topGO)
library(biomaRt)

#Read in data
mydata0<-read.table("OUfg_fullReproSet.selecAssocTrait.noErrors.txt", header=TRUE, sep="\t")
mydata<-mydata0[order(mydata0$protID),]

print("Converting mouse protein IDs to gene names...")
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org/")
all_genes<-getBM(attributes=c('ensembl_peptide_id','ensembl_gene_id','external_gene_name'), mart=ens_mus)
colnames(all_genes)<-c("prot_id","gene_id","gene_name")
all_genes_sorted<-all_genes[order(all_genes$prot_id),]
all_genes_uniq<-all_genes_sorted[which(!(duplicated(all_genes_sorted$prot_id))),]
all_genes_final<-all_genes_uniq[which(all_genes_uniq$prot_id != ""),]
all_genes_inDat<-all_genes_final[which(all_genes_final$prot_id %in% mydata$protID),]
print(dim(all_genes_inDat))
print(head(all_genes_inDat))
dat_hasGeneID<-mydata[which(mydata$protID %in% all_genes_inDat$prot_id),]
stopifnot(all.equal(dat_hasGeneID$protID, all_genes_inDat$prot_id))
mouse_gene<-unlist(sapply(dat_hasGeneID$protID, function(x) all_genes_inDat$gene_id[which(all_genes_inDat$prot_id==x)]))
stopifnot(length(mouse_gene)==nrow(dat_hasGeneID))
rownames(dat_hasGeneID)<-mouse_gene
print(paste("Number of genes lost in converting to gene ID:", nrow(mydata)-nrow(dat_hasGeneID)))

print("Running topGO analysis...")
subset<-rownames(dat_hasGeneID)[which(dat_hasGeneID$selec_result == "Selection is associated with the phenotype / trait")]
print(paste("Number of genes in focal subset:", length(subset)))
print(paste("Total number of genes in dataset:", nrow(dat_hasGeneID)))
total<-rownames(dat_hasGeneID)
inSubset<-c()
for(gene in total){
	if(gene %in% subset){
		inSubset<-c(inSubset, 1)
	} else{
		inSubset<-c(inSubset, 0)
	}
}
if(length(subset) != length(which(inSubset==1))){
	print("ERROR: length of 'subset' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
	return(iO==1)
}
myGOdata<-new("topGOdata", description="selAssocTraitVtotal", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("OUfg_fullReproSet.selecAssocTrait.GOenrichment.csv"), row.names=FALSE)
