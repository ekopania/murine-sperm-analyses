#PURPOSE: Get only selected (max targets) transcripts

mydata<-read.table("../repro_exons.protein_coding_only.exon_counts.txt", header=TRUE)
print(head(mydata))
selec_txpts<-scan("../selected-transcripts-targets.txt", what=character(), comment.char = "#")
print(head(selec_txpts))

newdata<-c()
for(i in 1:nrow(mydata)){
	if(mydata$ensembl_transcript_id[i] %in% selec_txpts){
		newdata<-rbind(newdata, mydata[i,])
	}
}

newdf<-as.data.frame(newdata)
print(dim(newdf))
write.table(newdf, file="../repro_exons.selected_txpts.exon_counts.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
