#PURPOSE: Plot median dN/dS for different groups of genes and compare (violin plots + Wilcoxon test)
#For dN/dS calculatd using PAML codeml M0 model

#library(ggplot2)

#Which busted run to extract data from
#Options: "" (SRV and MNMs), "noSRV_noMNM"", "yesSRV_noMNM" 
busted_run<-"noSRV_noMNM"

#Read in dN/dS
dnds<-read.csv(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet",busted_run,"csv", sep="."), header=TRUE, comment.char="#")
dnds$protID<-unlist(sapply(dnds$file, function(x) gsub("-mafft-cds.filter.json", "", x)))

#Read in gene sets for tissues/cell types
som<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

bd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BD.strict.txt", what=character())
bu<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BU.strict.txt", what=character())
cg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.CG.strict.txt", what=character())
dp<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.DP.strict.txt", what=character())
vp<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.VP.strict.txt", what=character())
sv<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.SV.strict.txt", what=character())

ts<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.testisSpecific.txt", what=character())

#Medians and Wilcoxon test
print(paste("Median dN/dS for all genes:", median(dnds$dn.ds, na.rm=TRUE)))
print(paste("Median rate of double nuc mut for all genes:", median(dnds$mnm2, na.rm=TRUE)))
print(paste("Median rate of triple nuc mut for all genes:", median(dnds$mnm3, na.rm=TRUE)))
print(paste("Number of genes in dataset:", nrow(dnds)))

groups<-c("som", "spg", "pre", "spc", "spd", "elo", "bd", "bu", "cg", "dp", "vp", "sv", "ts")
med_dnds<-c()
ns<-c()
wilcox_ps<-c()
med_mnm2<-c()
med_mnm3<-c()
for(i in groups){
	#Compare dN/dS
	print(paste("Median dN/dS for", i, "genes:", median(dnds$dn.ds[which(dnds$protID %in% get(i))], na.rm=TRUE)))
	print(paste("Number of", i, "genes in dataset:", length(which(dnds$protID %in% get(i)))))
	med_dnds<-c(med_dnds, median(dnds$dn.ds[which(dnds$protID %in% get(i))], na.rm=TRUE))
	ns<-c(ns, length(which(dnds$protID %in% get(i))))
	result<-wilcox.test(dnds$dn.ds, dnds$dn.ds[which(dnds$protID %in% get(i))])
	print(result)
	wilcox_ps<-c(wilcox_ps, result$p.value)

	#Repeat for MNMs
	print(paste("Median rate of double nuc mut for", i, "genes:", median(dnds$mnm2[which(dnds$protID %in% get(i))], na.rm=TRUE)))
	med_mnm2<-c(med_mnm2, median(dnds$mnm2[which(dnds$protID %in% get(i))], na.rm=TRUE))
	print(paste("Median rate of triple nuc mut for", i, "genes:", median(dnds$mnm3[which(dnds$protID %in% get(i))], na.rm=TRUE)))
        med_mnm3<-c(med_mnm3, median(dnds$mnm3[which(dnds$protID %in% get(i))], na.rm=TRUE))
}

print("Benjamini-Hochberg corrected p-values for Wilcoxon rank sum tests:")
wilcox_ps_adj<-p.adjust(wilcox_ps, method="BH")
names(wilcox_ps_adj)<-groups
print(wilcox_ps_adj)

final_df<-as.data.frame(rbind(med_mnm2, med_mnm3, med_dnds, ns, wilcox_ps, wilcox_ps_adj))
colnames(final_df)<-groups
rownames(final_df)<-c("mnm2", "mnm3", "dN/dS", "num_genes", "pval", "p.adj")
write.csv(final_df, file=paste("dnds_results_busted",busted_run,"csv", sep="."), quote=FALSE, row.names=TRUE)
