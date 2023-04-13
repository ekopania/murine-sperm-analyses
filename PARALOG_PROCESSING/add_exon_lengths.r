#PURPOSE: Add exon lengths to exon positive selection output

exons<-read.table("../SELEC_TESTS/pos_selec.pValCor.m1aVm2a.RTM_byExon.withProtID.txt", header=TRUE)
myinfo<-read.table("../repro_exons_full.txt", header=TRUE)
myinfo_exons<-myinfo[which(myinfo$ensembl_exon_id %in% exons$ensembl_exon_id),]
sorted_info<-myinfo_exons[order(myinfo_exons$external_gene_name),]
sorted_info_filt<-sorted_info[which(sorted_info$ensembl_peptide_id %in% exons$ensembl_peptide_id),]
sorted_info_filt_moresort<-sorted_info_filt[order(sorted_info_filt$ensembl_exon_id),]
sorted_info_filt_moresort<-sorted_info_filt[order(sorted_info_filt$external_gene_name, sorted_info_filt$ensembl_exon_id),]
print(head(sorted_info_filt_moresort))
print(head(exons))
stopifnot(all.equal(exons$ensembl_exon_id, sorted_info_filt_moresort$ensembl_exon_id))
stopifnot(all.equal(exons$ensembl_peptide_id, sorted_info_filt_moresort$ensembl_peptide_id))
stopifnot(all.equal(exons$external_gene_name, sorted_info_filt_moresort$external_gene_name))
lens<-sorted_info_filt_moresort$stops - sorted_info_filt_moresort$starts

newdata<-as.data.frame(cbind(exons, exon_lengths=lens))
print(head(newdata))

write.table(newdata, "../SELEC_TESTS/pos_selec.pValCor.m1aVm2a.RTM_byExon.withProtID.withLengths.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
