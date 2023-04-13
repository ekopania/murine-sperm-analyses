#PURPOSE: Make a table of gene names, gene IDs, transcript IDs, peptide IDs, and exon IDs for all reproductive genes we are interested in
library(biomaRt)

#Get gene info from Ensembl
#Host is important for keeping versions consistent!!!
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org/")
all_genes<-getBM(attributes=c('external_gene_name','ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id','ensembl_exon_id','chromosome_name','5_utr_start','5_utr_end','3_utr_start','3_utr_end'), mart=ens_mus)
print(dim(all_genes))
print(head(all_genes))

#Read in reproductive genes of interest
#myGenes<-scan("../repro_genes.txt", what=character())
myGenes<-scan("/uufs/chpc.utah.edu/common/HIPAA/proj_clarklab/kopania/temp_murinae_paralogs/repro_genes_full.txt")
print(myGenes)

#Extract reproductive genes from big table
repro_genes_table<-all_genes[which(all_genes$external_gene_name %in% myGenes),]
print(dim(repro_genes_table))
print(head(repro_genes_table))
print(tail(repro_genes_table))

#Save info from reproductive genes
write.table(repro_genes_table, file="../repro_genes_table.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
