#PURPOSE: Extract lists of tissue-enriched genes from Dean et al. 2009 (MBE)

library(biomaRt)

#Read in data
f1<-read.csv("deanEtal2009_mbe-08-0837-File001_msp094.csv", header=TRUE)
f2<-read.csv("deanEtal2009_mbe-08-0837-File002_msp094.csv", header=TRUE)

#Get genes in both datasets and get protein IDs
geneIDs<-union(f1$gene, f1$ENSMUSG)
print(paste("Length of gene ID vector:", length(geneIDs)))

#Get the tissue each gene shows maximum expression in
max_tissues<-c()
for(g in geneIDs){
	mt1<-f1$max_tissue[which(f1$gene==g)]

	#Make sure max tissues consistent across datasets
	if(g %in% f2$ENSMUSG){
		mt2<-f2$max_tissue[which(f2$ENSMUSG==g)]
		if(mt1==mt2){
			max_tissues<-c(max_tissues, mt1)
		} else{
			print(paste("Warning: max tissue not the same across datasets for", g))
			max_tissues<-c(max_tissues, NA)
		}
	} else{
		max_tissues<-c(max_tissues, mt1)
	}
}
names(max_tissues)<-geneIDs

#Load biomaRt
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','ensembl_peptide_id'), mart=ens_mus)
colnames(all_genes)<-c("protein_id", "gene_id")
print(head(all_genes))

#Loop through all tissues
#NOTE: In the paper there is "AP" (anterior prostate) but in the supp files there is no AP, but there is CG
#I suspect CG is coagulating gland and that they're using CG and AP interchangeably?
myTissues<-c("SV","CG","VP","DP","BU","BD")
for(t in myTissues){
	print(paste("Working on tissue:", t))
	#Get gene lists by tissue
	geneIDs_strict<-names(max_tissues)[which(max_tissues==t)]
	#Change gene ID to protein IDs
	protIDs_strict<-all_genes$protein_id[which(all_genes$gene_id %in% geneIDs_strict)]
	#Get unique protein IDs only
	protIDs_strict_unique<-unique(protIDs_strict)
	#Remove ""
	protIDs_strict_filter<-protIDs_strict_unique[which(protIDs_strict_unique != "")]
	print(paste("Total number of protein IDs for this tissue (strict):", length(protIDs_strict_filter)))
	#Write output
	write(protIDs_strict_filter, file=paste("prot_list",t,"strict.txt", sep="."), ncolumns=1)

	#Get lists of genes detected in tissues at all
	geneIDs_any<-union(f1$gene[which(f1[,t]>0)], f2$ENSMUSG[which(f2[,t]>0)])
	#Convert to protein IDs, unique, remove ""
	protIDs_any<-all_genes$protein_id[which(all_genes$gene_id %in% geneIDs_any)]
	protIDs_any_unique<-unique(protIDs_any)
	protIDs_any_filter<-protIDs_any_unique[which(protIDs_any_unique != "")]
	print(paste("Total number of protein IDs for this tissue (any):", length(protIDs_any_filter)))
        write(protIDs_any_filter, file=paste("prot_list",t,"any.txt", sep="."), ncolumns=1)
}

print("Done with 00_getDeanEtal_proteins.r")
