#PURPOSE: Get biotypes for all peptide IDs associated with reproductive genes; only keep protein coding genes

library(biomaRt)

mydata<-read.table("../repro_exons_full.with_exon_counts.txt", header=TRUE)

uniq_pid<-unique(sort(mydata$ensembl_peptide_id))
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="nov2020.archive.ensembl.org")
geneInfo<-getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'ensembl_exon_id', 'gene_biotype', 'transcript_biotype'), filters='ensembl_peptide_id', values=uniq_pid, mart=ens_mus)

repro_geneInfo<-geneInfo[which(geneInfo$ensembl_peptide_id %in% uniq_pid),]

keep<-c()
for(i in 1:nrow(repro_geneInfo)){
	if(repro_geneInfo$transcript_biotype[i] == "protein_coding"){
		keep<-c(keep, repro_geneInfo$ensembl_peptide_id[i])
	}
}

keep_data<-mydata[which(mydata$ensembl_peptide_id %in% keep),]

print(dim(keep_data))
print(head(keep_data))
write.table(keep_data, file="../repro_exons.protein_coding_only.exon_counts.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
