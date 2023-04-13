#PURPOSE: For a given list of genes, determine where in the filtering process it was removed from the dataset (if at all)

library(biomaRt)

args<-commandArgs(TRUE)
if(length(args)!=1){
        stop("Enter the name of a file containing a list of ensembl protein IDs as a command line argument")
}
myGenes<-args[1]

#Text file containing gene or list of genes to check, in ensemble peptide ID format (ENSMUSP000...)
#myGenes<-"prot_list.SV.any.txt"
query_set<-scan(myGenes, what=character())
#OR csv with genes to check in a column titled "ensPID"
#myGenes<-"../RERconverge/interesting_repro_genes.csv"
#allGenes<-read.csv(myGenes, header=TRUE)
#query_set<-allGenes$ensPID

print(paste("List of proteins:", myGenes))
print(paste("There are", length(query_set), "genes in the query set"))

#Append gene names to protein IDs
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl") #, host = "https://dec2021.archive.ensembl.org/")
all_genes<-getBM(attributes=c('ensembl_peptide_id','ensembl_gene_id','external_gene_name'), mart=ens_mus)
colnames(all_genes)<-c("prot_id","gene_id","gene_name")
query_gID<-unlist(sapply(query_set, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))
query_gName<-unlist(sapply(query_set, function(x) all_genes$gene_name[which(all_genes$prot_id == x)]))
query_df<-as.data.frame(cbind(query_set, query_gID, query_gName))
print(query_df)

#Read in files used for checking filtering steps
probe_set<-scan("/mnt/beegfs/ek112884/murinae/EXOME_TARGETS/prot_list.all_targets.sort_uniq.txt", what=character())
orthos<-read.table("/mnt/beegfs/gt156213e/murinae-seq/02-Annotation/mouse-rat-orths-ens99.txt", header=TRUE, sep="\t")
print(head(orthos))
max_targets<-read.table("/mnt/beegfs/gt156213e/murinae-seq/02-Annotation/max-target-transcripts/selected-transcripts-targets.txt", header=TRUE, sep="\t")
pre_trim<-scan("prot_list.full-coding-mafft-f150.txt", what=character())
final_set<-scan("prot_list.full-coding-trimmed-f175-seq20-site50.txt",, what=character())

#Loop through all genes and save reason each gene was filtered out (if at all)
filt<-c()
for(g in query_set){
	#print(g)
	#Get row number from ortholog table for this protein ID
	thisRow<-which(orthos$Protein.stable.ID==g)
	#print(orthos[thisRow,])
	#print(length(thisRow))
	#Check if in probes
	if(!(g %in% probe_set)){
		filt<-c(filt, "no_targets")
	} else if(length(thisRow) == 0){
		filt<-c(filt, "not_in_ortho_file")
	} else if(orthos$Rat.gene.stable.ID[thisRow] == ""){
		filt<-c(filt, "no_rat_ortho")
	} else if(orthos$Rat.homology.type[thisRow] != "ortholog_one2one"){
		filt<-c(filt, "ortho_not_one2one")
	} else if(orthos$Rat.orthology.confidence..0.low..1.high.[thisRow] != 1){
		filt<-c(filt, "low_ortho_conf")
	} else if(is.na(orthos$dS.with.Rat[thisRow])){ #Need to start with this or the > will cause an error if NA
		filt<-c(filt, "high_dS")
	} else if(orthos$dS.with.Rat[thisRow] > 0.5){
		filt<-c(filt, "high_dS")
	} else if(!(g %in% max_targets$Protein.stable.ID)){
		filt<-c(filt, "not_longest_targets")
	} else if( (g %in% pre_trim) && (!(g %in% final_set)) ){
		filt<-c(filt, "trimming")
	} else{
		filt<-c(filt, NA)
	}
}

names(filt)<-query_set

print("BY TRANSCRIPT:")
print(filt)
print(paste("Number not in target probes:", length(which(filt == "no_targets"))))
print(paste("Number not in ortholog file:", length(which(filt == "not_in_ortho_file"))))
print(paste("Number with no rat ortholog:", length(which(filt == "no_rat_ortho"))))
print(paste("Number with no 1:1 rat ortholog:", length(which(filt == "ortho_not_one2one"))))
print(paste("Number with no high confidence rat ortholog:", length(which(filt == "low_ortho_conf"))))
print(paste("Number with high dS compared to rat:", length(which(filt == "high_dS"))))
print(paste("Number not maximum targets within gene:", length(which(filt == "not_longest_targets"))))
print(paste("Number filtered out at trimming:", length(which(filt == "trimming"))))

print(paste("Number overlapping with final dataset:", length(intersect(query_set, final_set))))
print(paste("Number not filtered out:", length(which(is.na(filt)))))
print("^ these last two should be the same")

print("BY GENE")
filt_df<-as.data.frame(cbind(query_df, filt))
print(head(filt_df))
filt_df_sort<-filt_df[order(filt_df$query_gID),]
print(head(filt_df_sort))
filt_longest<-filt_df_sort[union(which(is.na(filt_df_sort$filt)), which(!(filt_df_sort$filt == "not_longest_targets"))),]
print(filt_longest)
filt_unique<-filt_longest[!(duplicated(filt_longest$query_gID)),]
print(filt_unique)

print(paste("Number not in target probes:", length(which(filt_unique$filt == "no_targets"))))
print(paste("Number not in ortholog file:", length(which(filt_unique$filt == "not_in_ortho_file"))))
print(paste("Number with no rat ortholog:", length(which(filt_unique$filt == "no_rat_ortho"))))
print(paste("Number with no 1:1 rat ortholog:", length(which(filt_unique$filt == "ortho_not_one2one"))))
print(paste("Number with no high confidence rat ortholog:", length(which(filt_unique$filt == "low_ortho_conf"))))
print(paste("Number with high dS compared to rat:", length(which(filt_unique$filt == "high_dS"))))
#print(paste("Number not maximum targets within gene:", length(which(filt == "not_longest_targets"))))
print(paste("Number filtered out at trimming:", length(which(filt_unique$filt == "trimming"))))

print(paste("Number overlapping with final dataset:", length(intersect(filt_unique$query_set, final_set))))
print(paste("Number not filtered out:", length(which(is.na(filt_unique$filt)))))
print("^ these last two should be the same")

warnings()
