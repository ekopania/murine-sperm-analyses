#PURPOSE: Plot proportion of genes under positive selection based on busted in HyPhy for spermatogenesis cell types
#Also plot propotion of genes that are testis specific 

library(ggplot2)
library(gridExtra)

#Which busted run to extract data from
#Options: "" (SRV and MNMs), "noSRV_noMNM"", "yesSRV_noMNM"
busted_run<-"noSRV_noMNM"

#Read in and merge dN/dS
posSelec<-read.csv(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet", busted_run,"csv", sep="."), header=TRUE, comment.char="#")
print(head(posSelec))
print(dim(posSelec))

posSelec_filtered<-posSelec
posSelec_filtered$file<-unlist(sapply(posSelec_filtered$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
posSelec_filtered$p.adj<-p.adjust(posSelec_filtered$pval, method="BH")
pos_selec_genes<-posSelec_filtered$file[which(posSelec_filtered$p.adj < 0.05)]
all_selec_genes<-posSelec_filtered$file
print(paste("Number of genes in analysis:", length(all_selec_genes)))
print(paste("Number of genes with evidence for positive selection:", length(pos_selec_genes)))
prop_selec_all<-length(pos_selec_genes) / length(all_selec_genes)
print(prop_selec_all)

#Read in gene sets for tissues/cell types
som<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

ts<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.testisSpecific.txt", what=character())

#Get some stats for all genes
prop_ts_all<-length(intersect(ts, all_selec_genes)) / length(all_selec_genes)
prop_ts_selec_all<-length(intersect(ts, pos_selec_genes))/length(intersect(ts, all_selec_genes))
print(paste(length(intersect(ts, pos_selec_genes)), "testis-specific are under positive selection out of", length(intersect(ts, all_selec_genes))))

#Get proportions and plot
groups<-c("som", "spg", "pre", "spc", "spd", "elo", "ts")
ns<-c()
prop_selec<-c()
num_selec<-c()
prop_test_ps<-c()
prop_ts<-c()
num_ts<-c()
prop_ts_selecSet<-c()
prop_ts_selecSet_ps<-c()
prop_selec_ts<-c()
prop_selec_ts_test_ps<-c()
for(i in groups){
	print(paste("Working on", i))
	#Positive selection
	test_set<-posSelec_filtered[which(posSelec_filtered$file %in% get(i)),]
	print(dim(test_set))
	ns<-c(ns, nrow(test_set))
	#print(test_set)
	print(paste("Number of genes in cell type with evidence for positive selection:", length(which(test_set$p.adj < 0.05))))
	prop<-length(which(test_set$p.adj < 0.05)) / nrow(test_set)
	prop_selec<-c(prop_selec, prop)
	num_selec<-c(num_selec, length(which(test_set$p.adj < 0.05)))
	result<-prop.test(length(which(test_set$p.adj < 0.05)), nrow(test_set), prop_selec_all)
	prop_test_ps<-c(prop_test_ps, result$p.value)
	#Testis-specific - all
	print(paste("Number of genes in cell type that are testis-specific:", length(intersect(get(i), ts))))
	pts<-length(intersect(get(i), ts)) / length(get(i))
	prop_ts<-c(prop_ts, pts)
	num_ts<-c(num_ts, length(intersect(get(i), ts)))
	#Testis-specific - in selection test gene set ONLY
	ts_selec_genes<-Reduce(intersect, list(get(i), ts, all_selec_genes))
	print(paste("Number of genes in cell type that are testis-specific and included in selection test dataset:", length(ts_selec_genes)))
	prop_ts_selecSet<-c(prop_ts_selecSet, length(ts_selec_genes)/nrow(test_set))
	result<-prop.test(length(ts_selec_genes), nrow(test_set), prop_ts_all)
	prop_ts_selecSet_ps<-c(prop_ts_selecSet_ps, result$p.value)
	#Positive selection and testis-specific
	ps_ts<-intersect(test_set$file[which(test_set$p.adj < 0.05)], ts)
	print(paste("Number of testis-specific genes in cell type with evidence for positive selection:", length(ps_ts)))
	if(length(ts_selec_genes) > 0){
		prop_selec_ts<-c(prop_selec_ts, length(ps_ts)/length(ts_selec_genes))
		result<-prop.test(length(ps_ts), length(ts_selec_genes), prop_ts_selec_all)
		prop_selec_ts_test_ps<-c(prop_selec_ts_test_ps, result$p.value)
	} else{
		prop_selec_ts<-c(prop_selec_ts, NA)
		prop_selec_ts_test_ps<-c(prop_selec_ts_test_ps, NA)
	}
}

print("Proportion of genes with evidence for positive selection:")
names(prop_selec)<-groups
print(prop_selec)

print("Proportion test results, compared to genome-wide positive selection:")
#names(prop_test_ps)<-groups
#print(prop_test_ps)
prop_test_ps_adj<-p.adjust(prop_test_ps, method="BH")
names(prop_test_ps_adj)<-groups
print(prop_test_ps_adj)

#print("Pairwise proportion tests for positive selection:")
#names(num_selec)<-c("som","spg","pre","spc","spd","elo","ts")
#pairwise.prop.test(num_selec[1:6], c(length(som), length(spg), length(pre), length(spc), length(spd), length(elo), length(ts)), p.adjust.method="fdr")

print("Proportion of genes that are testis-specific:")
names(prop_ts)<-groups
print(prop_ts)

#print("Pairwise proportion tests for testis-specificity:")
#names(num_ts)<-c("som","spg","pre","spc","spd","elo", "ts")
#pairwise.prop.test(num_ts[1:6], c(length(som), length(spg), length(pre), length(spc), length(spd), length(elo)), p.adjust.method="fdr")

print("Proportions of genes that are testis-specific - selection test dataset ONLY:")
names(prop_ts_selecSet)<-groups
print(prop_ts_selecSet)
print("Proportion test results for testis-specific and in selection dataset:")
prop_ts_selecSet_ps_adj<-p.adjust(prop_ts_selecSet_ps, method="BH")
names(prop_ts_selecSet_ps_adj)<-groups
print(prop_ts_selecSet_ps_adj)

print("Proportions of testis-specific genes that are under selection:")
names(prop_selec_ts)<-groups
print(prop_selec_ts)
print("Proportion test results for testis-specific genes under selection:")
prop_selec_ts_test_ps_adj<-p.adjust(prop_selec_ts_test_ps, method="BH")
names(prop_selec_ts_test_ps_adj)<-groups
print(prop_selec_ts_test_ps_adj)

#myDF0<-as.data.frame(cbind(prop.pos.selec=prop_selec, prop.testis.specific=prop_ts, dataset=rep("sc",length(prop_selec))))
myDF0<-as.data.frame(cbind(prop.pos.selec=prop_selec, prop.testis.specific=prop_ts_selecSet, dataset=rep("sc",length(prop_selec))))
myDF<-myDF0[which(rownames(myDF0) != "ts"),]
rownames(myDF)<-c("somatic","spermatogonia","pre-leptotene","spermatocytes","spermatids","elongating")
print(head(myDF))

print("Plotting...")
p<-ggplot(myDF, aes(x=rownames(myDF), y=as.numeric(as.character(prop.pos.selec)), group=dataset)) + geom_line(size=2) + geom_point(size=4)
p<-p + labs(title="Proportion of genes under positive selection", x="Cell Type", y="Proportion of genes")
p<-p + scale_x_discrete(limits=c("spermatogonia","pre-leptotene","spermatocytes","spermatids","elongating")) #+ ylim(0,0.015)
p<-p + theme_minimal() + theme(axis.text.y=element_text(size=18), axis.title=element_text(size=21))
p<-p + geom_hline(yintercept=length(pos_selec_genes)/length(all_selec_genes), linetype="dashed")

p_ts<-ggplot(myDF, aes(x=rownames(myDF), y=as.numeric(as.character(prop.testis.specific)), group=dataset)) + geom_line(size=2) + geom_point(size=4)
p_ts<-p_ts + labs(title="Proportion of genes that are testis-specific", x="Cell Type", y="Proportion of genes")
p_ts<-p_ts + scale_x_discrete(limits=c("spermatogonia","pre-leptotene","spermatocytes","spermatids","elongating"))
p_ts<-p_ts + theme_minimal() + theme(axis.text.y=element_text(size=18), axis.title=element_text(size=21))
p_ts<-p_ts + geom_hline(yintercept=prop_ts_all, linetype="dashed")

plots.list<-list(p,p_ts)
plots<-marrangeGrob(plots.list,nrow=2,ncol=1)
ggsave(paste("prop_pos_selec.busted",busted_run,"green_sc.pdf", sep="."), plots, width = 12, height = 7, units = "in")
dev.off()

final_df<-as.data.frame(rbind(ns, prop_selec, prop_test_ps, prop_test_ps_adj, prop_ts_selecSet, prop_ts_selecSet_ps, prop_ts_selecSet_ps_adj, prop_selec_ts, prop_selec_ts_test_ps, prop_selec_ts_test_ps_adj))
colnames(final_df)<-groups
rownames(final_df)<-c("num_genes", "proportion under pos selec", "pos selec pval", "pos selec BH pval", "proportion testis-specific", "testis-specific pval", "testis-specific BH pval", "proportion of testis-specific genes under pos selec", "testis-specific pos selec pval", "testis-specific pos selec BH pval")
write.csv(final_df, paste("selection_results.busted",busted_run,"cell_types.csv", sep="."), row.names=TRUE, quote=FALSE)

