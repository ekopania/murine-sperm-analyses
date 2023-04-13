#PURPOSE: Make a bedfile of exons in all reproductive genes we are interested in

library(biomaRt)

#Read in repro genes table
myData<-read.table("../repro_genes_table.txt", header=TRUE)

#Get gene info from Ensembl
#Host is important for keeping versions consistent!!!
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="nov2020.archive.ensembl.org")
all_genes<-getBM(attributes=c('chromosome_name','exon_chrom_start','exon_chrom_end','ensembl_exon_id','strand','phase'), mart=ens_mus)
print(dim(all_genes))
print(head(all_genes))

#Extract reproductive genes from big table
repro_genes_bed<-all_genes[which(all_genes$ensembl_exon_id %in% myData$ensembl_exon_id),]
print(dim(repro_genes_bed))
print(head(repro_genes_bed))
print(tail(repro_genes_bed))

#Make it zero based, non-inclusive (for bedfile)
repro_genes_bed$exon_chrom_start<-repro_genes_bed$exon_chrom_start - 1
print(dim(repro_genes_bed))
print(head(repro_genes_bed))
print(tail(repro_genes_bed))

#Save info from reproductive genes
write.table(repro_genes_bed, file="repro_exons.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
