#PURPOSE: Plot proportion of genes under positive selection based on busted in HyPHy for reproductive tissues
#Also plot propotion of genes that are testis specific 

library(ggplot2)
library(gridExtra)

#Which busted run to extract data from
#Options: "" (SRV and MNMs), "noSRV_noMNM"", "yesSRV_noMNM"
busted_run<-"noSRV_noMNM"

#Read in and merge dN/dS
posSelec<-read.csv(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet",busted_run,"csv", sep="."), header=TRUE, comment.char="#")

posSelec_filtered<-posSelec
posSelec_filtered$file<-unlist(sapply(posSelec_filtered$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
posSelec_filtered$p.adj<-p.adjust(posSelec_filtered$pval, method="BH")
pos_selec_genes<-posSelec_filtered$file[which(posSelec_filtered$p.adj < 0.05)]
all_selec_genes<-posSelec_filtered$file
print(paste("Number of genes in analysis:", length(all_selec_genes)))
print(paste("Number of genes with evidence for positive selection:", length(pos_selec_genes)))
prop_selec_all<-length(pos_selec_genes) / length(all_selec_genes)

#Read in gene sets for tissues/cell types
bd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BD.any.txt", what=character())
bu<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BU.any.txt", what=character())
cg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.CG.any.txt", what=character())
dp<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.DP.any.txt", what=character())
vp<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.VP.any.txt", what=character())
sv<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.SV.any.txt", what=character())

#Get proportions and plot
groups<-c("bd", "bu", "cg", "dp", "vp", "sv")
ns<-c()
prop_selec<-c()
prop_test_ps<-c()
for(i in groups){
	print(paste("Working on", i))
	test_set<-posSelec_filtered[which(posSelec_filtered$file %in% get(i)),]
	print(dim(test_set))
	ns<-c(ns, nrow(test_set))
	print(paste("Number of genes in cell type with evidence for positive selection:", length(which(test_set$p.adj < 0.05))))
	prop<-length(which(test_set$p.adj < 0.05)) / nrow(test_set)
	prop_selec<-c(prop_selec, prop)
	result<-prop.test(length(which(test_set$p.adj < 0.05)), nrow(test_set), prop_selec_all)
        prop_test_ps<-c(prop_test_ps, result$p.value)
}

print("Proportion of genes with evidence for positive selection:")
names(prop_selec)<-groups
print(prop_selec)

print("Proportion test results, compared to genome-wide positive selection:")
names(prop_test_ps)<-groups
print(prop_test_ps)
prop_test_ps_adj<-p.adjust(prop_test_ps, method="BH")
names(prop_test_ps_adj)<-groups
print(prop_test_ps_adj)

myDF<-as.data.frame(cbind(prop.pos.selec=prop_selec, groups=rep("tissue", length(prop_selec))))
rownames(myDF)<-c("bulbourethral_diverticulum","bulbourethral_gland","coagulating_gland","dorsal_prostate","ventral_prostate","seminal_vesicle")
print(head(myDF))

print("Plotting...")
p<-ggplot(myDF, aes(x=rownames(myDF), y=as.numeric(as.character(prop.pos.selec)), group=groups)) + geom_line(size=2) + geom_point(size=4)
p<-p + labs(title="Proportion of genes under positive selection", x="", y="Proportion of genes")
p<-p + scale_x_discrete(limits=c("bulbourethral_diverticulum","bulbourethral_gland","coagulating_gland","dorsal_prostate","ventral_prostate","seminal_vesicle")) #+ ylim(0,0.015)
p<-p + theme_minimal() + theme(axis.text.y=element_text(size=18), axis.title=element_text(size=21))
p<-p + geom_hline(yintercept=length(pos_selec_genes)/length(all_selec_genes), linetype="dashed")

pdf(paste("prop_pos_selec.busted",busted_run,"dean_tissues.pdf", sep="."))
print(p)
dev.off()

final_df<-as.data.frame(rbind(ns, prop_selec, prop_test_ps, prop_test_ps_adj))
colnames(final_df)<-groups
rownames(final_df)<-c("num_genes", "proportion under pos selec", "pos selec pval", "pos selec BH pval")
write.csv(final_df, paste("selection_results.busted",busted_run,"tissues.csv", sep="."), row.names=TRUE, quote=FALSE)
