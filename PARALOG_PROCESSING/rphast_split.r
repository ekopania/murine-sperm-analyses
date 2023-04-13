#PURPOSE: Use rphast sub.msa to split alignments by exon

library(rphast)

#I think I'll need to do this manually for each gene since exon coordinates will be different for each one

#Ovgp1
#ENSMUSP00000000573; ENSMUSP00000000573; start 292; stop 1297
#Svs1

#Svs2
#ENSMUSP00000042389; ENSMUSE00000244955; start 1; stop 1051

#Svs3a
#ENSMUSP00000104995; ENSMUSE00000982347
#ENSMUSP00000017147; ENSMUSE00000982347

#Svs4
#ENSMUSP00000017142; ENSMUSE00000171642;
#ENSMUSP00000017142; ENSMUSE00000171641; 
#ENSMUSP00000017142; ENSMUSE00000406734

#Svs5
#ENSMUSP00000017148; ENSMUSE00000171645
#ENSMUSP00000119963; ENSMUSE00000171645

#Svs6
#ENSMUSP00000017144; ENSMUSE00000244879
#ENSMUSP00000017144; ENSMUSE00000161934
#ENSMUSP00000017144; ENSMUSE00000161935
#ENSMUSP00000017144; ENSMUSE00000244864

protein<-"ENSMUSP00000042389"
exon<-"ENSMUSE00000244955"
start<-1
stop<-1051

#Read in fasta alignment
in_msa<-read.msa(paste0("../MAFFT-f175/nt/",protein,"-mafft-cds.fa"), format="FASTA")

#Get start and end columns for exon
mm10<-sub.msa(in_msa, seqs="mm10", keep=TRUE)
myseq<-as.character(mm10$seq)
myseq_split<-strsplit(myseq, "")[[1]]
myseq_sep<-c()
for(i in myseq_split){
	myseq_sep<-c(myseq_sep, i)
}
gap_index<-which(myseq_sep != "-")

my_sub<-sub.msa(x=in_msa, start.col=start, end.col=stop, refseq="mm10")

write.msa(my_sub, file=paste(protein,exon,"submsa.fa", sep="_"), format="FASTA")
