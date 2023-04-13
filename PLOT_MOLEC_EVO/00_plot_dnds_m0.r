#PURPOSE: Plot median dN/dS for different groups of genes and compare (violin plots + Wilcoxon test)
#For dN/dS calculatd using PAML codeml M0 model

library(ggplot2)

#Read in dN/dS
dnds<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/RTM_SET_pared_M0/dnds.RTM_SET_pared_M0.txt", header=TRUE)

#Read in gene sets for tissues/cell types
som<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

bd<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BD.any.txt", what=character())
bu<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BU.any.txt", what=character())
cg<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.CG.any.txt", what=character())
dp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.DP.any.txt", what=character())
vp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.VP.any.txt", what=character())
sv<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.SV.any.txt", what=character())

ts<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_testisSpecific.txt", what=character())

#Medians and Wilcoxon test
print(paste("Median dN/dS for all genes:", median(dnds$dN.dS, na.rm=TRUE)))
print(paste("Number of genes in dataset:", nrow(dnds)))

groups<-c("som", "spg", "pre", "spc", "spd", "elo", "bd", "bu", "cg", "dp", "vp", "sv", "ts")
wilcox_ps<-c()
for(i in groups){
	print(paste("Median dN/dS for", i, "genes:", median(dnds$dN.dS[which(dnds$protID %in% get(i))], na.rm=TRUE)))
	print(paste("Number of", i, "genes in dataset:", length(which(dnds$protID %in% get(i)))))
	result<-wilcox.test(dnds$dN.dS, dnds$dN.dS[which(dnds$protID %in% get(i))])
	print(result)
	wilcox_ps<-c(wilcox_ps, result$p.value)
}

print("Benjamini-Hochberg corrected p-values for Wilcoxon rank sum tests:")
wilcox_ps_adj<-p.adjust(wilcox_ps, method="BH")
names(wilcox_ps_adj)<-groups
print(wilcox_ps_adj)
