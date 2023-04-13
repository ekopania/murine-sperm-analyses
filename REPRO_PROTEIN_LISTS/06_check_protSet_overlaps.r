#PURPOSE: Check if given gene or set of genes was present in exon capture probe set

#Text file containing gene or list of genes to check, in ensemble peptide ID format (ENSMUSP000...)
#myGenes<-"prot_list.SV.strict.txt"
myGenes<-"../RERconverge/interesting_repro_genes.csv"
allGenes<-read.csv(myGenes, header=TRUE)
#target_dataset<-"prot_list.reproductive-testes-mass-coding-trimmed-f0-seq20-site50.txt"
#target_dataset<-"prot_list.full-coding-trimmed-f175-seq20-site50.txt"
target_dataset<-"prot_list.full-coding-mafft-f150.txt"
print(paste("List of proteins:", myGenes))
print(paste("List of target proteins:", target_dataset))

#query_set<-scan(myGenes, what=character())
query_set<-allGenes$ensPID
probe_set<-scan("/mnt/beegfs/ek112884/murinae/EXOME_TARGETS/prot_list.all_targets.sort_uniq.txt", what=character())
target_set<-scan(target_dataset, what=character())

print("Overlap with exon capture probes:")
print(length(intersect(query_set, probe_set)))
print(intersect(query_set, probe_set))

print("In query but NOT in probe set:")
print(length(which(!(query_set %in% probe_set))))
print(query_set[which(!(query_set %in% probe_set))])

print("Overlap with target dataset:")
print(length(intersect(query_set, target_set)))
print(intersect(query_set, target_set))
