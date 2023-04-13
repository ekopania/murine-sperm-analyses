#PURPOSE: Append gene symbols and peptide IDs to pos selection results for runs done with exons, instead of proteins

mydata<-read.table("../SELEC_TESTS/pos_selec.pValCor.m1aVm2a.RTM_byExon.txt", header=TRUE)
repro<-read.table("../repro_genes_table.txt", header=TRUE)
gene_txp_pep<-repro[,c("external_gene_name","ensembl_transcript_id", "ensembl_peptide_id","ensembl_exon_id")]

gp_append<-c()
for(i in 1:nrow(mydata)){
	this_exon<-mydata$protID[i]
	this_gp<-which(gene_txp_pep$ensembl_exon_id == this_exon)
	gp_append<-rbind(gp_append, gene_txp_pep[this_gp,])
}
print(dim(gp_append))

selec_txp<-scan("../selected-transcripts-targets.txt", what=character(), comment="#")
gp_append_select<-gp_append[which(gp_append$ensembl_transcript_id %in% selec_txp),]
print(dim(gp_append_select))
print(head(gp_append_select))
print(head(mydata))

stopifnot(all.equal(mydata$protID, gp_append_select$ensembl_exon_id))

newdata<-as.data.frame(cbind(gp_append_select, mydata[,2:6]))
sortdata<-newdata[order(newdata$external_gene_name),]
print(head(sortdata))

write.table(sortdata, file="../SELEC_TESTS/pos_selec.pValCor.m1aVm2a.RTM_byExon.withProtID.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


q()

dim(mydata)
[1] 123   6
tail(gp_append_select)
       external_gene_name ensembl_transcript_id ensembl_peptide_id
496117            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
496118            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
496119            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
496120            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
496121            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
496122            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
          ensembl_exon_id
496117 ENSMUSE00000158118
496118 ENSMUSE00000158118
496119 ENSMUSE00000158118
496120 ENSMUSE00000158118
496121 ENSMUSE00000158118
496122 ENSMUSE00000158118
tail(gp_append)
       external_gene_name ensembl_transcript_id ensembl_peptide_id
496120            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
507120            Slco6c1    ENSMUST00000189547 ENSMUSP00000140791
496121            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
507121            Slco6c1    ENSMUST00000189547 ENSMUSP00000140791
496122            Slco6c1    ENSMUST00000027569 ENSMUSP00000027569
507122            Slco6c1    ENSMUST00000189547 ENSMUSP00000140791
          ensembl_exon_id
496120 ENSMUSE00000158118
507120 ENSMUSE00000158118
496121 ENSMUSE00000158118
507121 ENSMUSE00000158118
496122 ENSMUSE00000158118
507122 ENSMUSE00000158118
>

