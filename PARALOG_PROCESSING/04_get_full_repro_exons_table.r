#PURPOSE: Merge annotation and bedfile info to generate a table with all exon information and coordinates

#Read in data
mydata<-read.table("../repro_genes_table.txt", header=TRUE)
print(dim(mydata))
print(head(mydata))

mybed<-read.table("../repro_exons.bed", header=FALSE)
print(dim(mybed))
print(head(mybed))

#targets<-read.table("/mnt/beegfs/ek112884/murinae/EXOME_TARGETS/mm9_to_mm10_targeted_regions.bed", header=FALSE)
#print(head(targets))

#Add start pos, stop pos, and strand from bedfile to annotaton table
#Remove exons that weren't' in capture set
starts<-c()
stops<-c()
strands<-c()
non_proco<-c()
no_bed<-c()
#no_targets<-c()

for(i in 1:nrow(mydata)){
	e<-mydata$ensembl_exon_id[i]
	if(is.na(mydata$ensembl_peptide_id[i])){
		print(paste("Non-protein-coding:", mydata$ensembl_peptide_id[i], e))
		non_proco<-c(non_proco, i)
	} else if(e %in% mybed$V4){
		strand<-as.character(mybed$V5[which(mybed$V4 == e)])
		if(!(is.na(mydata$X5_utr_end[i]))){
			if(strand == "+"){
				sta<-mydata$X5_utr_end[i] + 1 #Forward direction (strand + or 1)
			} else{
				sto<-mydata$X5_utr_start[i] - 1 #Reverse direction (strand - or -1); sto so that sto > sta always
			}
		} else{
			if(strand == "+"){
				sta<-mybed$V2[which(mybed$V4 == e)] + 1 #Forward; +1 to deal with 0 vs 1 base shenanigans
			} else{
				sto<-mybed$V3[which(mybed$V4 == e)] #Reverse
			}
		}
		if(!(is.na(mydata$X3_utr_start[i]))){
			if(strand == "+"){
                                sto<-mydata$X3_utr_start[i] - 1 #Forward
			} else{
				sta<-mydata$X3_utr_end[i] + 1 #Reverse
			}
		} else{
			if(strand == "+"){
				sto<-mybed$V3[which(mybed$V4 == e)] #Forward
			} else{
				sta<-mybed$V2[which(mybed$V4 == e)] + 1 #Reverse; +1 to deal with 0 vs 1 base shenanigans
			}
		}
		#if(sta %in% targets$V2){
			#if(targets$V3[which(targets$V2 == sta)] == sto){
				starts<-c(starts, sta)
				stops<-c(stops, sto)
				strands<-c(strands, strand)
			#} else{
			#	print(paste("Not in targets:", e, sta, sto))
	                #        no_targets<-c(no_targets, e)
			#}
		#} else{
		#	print(paste("Not in targets:", e, sta, sto))
		#	no_targets<-c(no_targets, e)
		#}
	} else{
		print(paste("Not in bedfile:", e))
		no_bed<-c(no_bed, i)
	}
	#print(length(starts))
	#print(length(stops))
	#print(length(strands))
	#print(length(no_bed))
	if(length(starts) != length(stops)){
		print("Problem exon:")
		print(e)
		print(length(starts))
		print(length(stops))
	}
}

print(length(starts))
print(length(stops))

#Remove exons that aren't in bed file 
filtdata1<-mydata[-(c(no_bed, non_proco)),]
#Remove exons that have no coding sequence
newdata<-as.data.frame(cbind(filtdata1, starts, stops, strands))
print(head(newdata))
print(dim(newdata))

diff<-newdata$stops - newdata$starts
print("Any starts greater than stops?")
print(newdata[which(diff < -1),])
filtdata2<-newdata[which(diff > 1),]
print(head(filtdata2))
print(dim(filtdata2))

#Write table
write.table(filtdata2, file="../repro_exons_full.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
