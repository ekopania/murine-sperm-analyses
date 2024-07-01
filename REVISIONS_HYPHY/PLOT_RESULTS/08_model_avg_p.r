#PURPOSE: Calculate model-averaged p-values from Lucaci et al. 2023 MBE

#Which tissue gene sets to use
#Options: any (detected in tissue) or strict (detected at highest level in tissue)
tissue_stringency<-"strict"

#Read in data from all BUSTED runs
srv_dtMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.yesSRV_dtMNM.csv", header=TRUE, comment.char="#")
srv_noMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.yesSRV_noMNM.csv", header=TRUE, comment.char="#")
noSRV_noMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.noSRV_noMNM.csv", header=TRUE, comment.char="#")

srv_dtMNM$protID<-unlist(sapply(srv_dtMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
srv_noMNM$protID<-unlist(sapply(srv_noMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
noSRV_noMNM$protID<-unlist(sapply(noSRV_noMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))

#Get loci that converged in all models
srv_dtMNM_conv<-srv_dtMNM$protID[which(srv_dtMNM$lrt >= 0)]
srv_noMNM_conv<-srv_noMNM$protID[which(srv_noMNM$lrt >= 0)]
noSRV_noMNM_conv<-noSRV_noMNM$protID[which(noSRV_noMNM$lrt >= 0)]
keep<-Reduce(intersect, list(srv_dtMNM_conv, srv_noMNM_conv, noSRV_noMNM_conv))

srv_dtMNM_keep<-srv_dtMNM[which(srv_dtMNM$protID %in% keep),]
srv_noMNM_keep<-srv_noMNM[which(srv_noMNM$protID %in% keep),]
noSRV_noMNM_keep<-noSRV_noMNM[which(noSRV_noMNM$protID %in% keep),]

srv_dtMNM_sort<-srv_dtMNM_keep[order(srv_dtMNM_keep$protID),]
srv_noMNM_sort<-srv_noMNM_keep[order(srv_noMNM_keep$protID),]
noSRV_noMNM_sort<-noSRV_noMNM_keep[order(noSRV_noMNM_keep$protID),]

stopifnot(all.equal(srv_dtMNM_sort$protID, srv_noMNM_sort$protID))
stopifnot(all.equal(srv_dtMNM_sort$protID, noSRV_noMNM_sort$protID))

#Make a new df with AIC and p-val for each model
sigdf<-as.data.frame(cbind(srv_dtMNM_sort$aic, srv_dtMNM_sort$pval, srv_noMNM_sort$aic, srv_noMNM_sort$pval, noSRV_noMNM_sort$aic, noSRV_noMNM_sort$pval))
colnames(sigdf)<-c("srv_dtMNM_aic", "srv_dtMNM_pval", "srv_noMNM_aic", "srv_noMNM_pval", "noSRV_noMNM_aic", "noSRV_noMNM_pval")

#Loop through each gene and calculate model-averaged p-val
ma_ps<-c()
for(i in 1:nrow(sigdf)){
#for(i in 1:10){
	#Determine best AIC of the models
	aic_best<-min(c(sigdf$srv_dtMNM_aic[i], sigdf$srv_noMNM_aic[i], sigdf$noSRV_noMNM_aic[i]))
	#Get Akaike weight for each model
	w_srv_dtMNM<-exp((aic_best - sigdf$srv_dtMNM_aic[i])/2)
	w_srv_noMNM<-exp((aic_best - sigdf$srv_noMNM_aic[i])/2)
	w_noSRV_noMNM<-exp((aic_best - sigdf$noSRV_noMNM_aic[i])/2)
	#Get weighted p-val for each model
	wp_srv_dtMNM<-w_srv_dtMNM * sigdf$srv_dtMNM_pval[i]
	wp_srv_noMNM<-w_srv_noMNM * sigdf$srv_noMNM_pval[i]
	wp_noSRV_noMNM<-w_noSRV_noMNM * sigdf$noSRV_noMNM_pval[i]
	#Sum weighted p-vals to get model-averaged p-val
	pma<-sum(wp_srv_dtMNM, wp_srv_noMNM, wp_noSRV_noMNM)
	ma_ps<-c(ma_ps, pma)
	#Print everything to check calculations
	#print(sigdf[i,])
	#print(aic_best)
	#print(w_srv_dtMNM)
	#print(w_srv_noMNM)
	#print(w_noSRV_noMNM)
	#print(wp_srv_dtMNM)
	#print(wp_srv_noMNM)
	#print(wp_noSRV_noMNM)
	#print(pma)
}

#Correct for multiple tests; add protein IDs
ma_p_adj<-p.adjust(ma_ps, method="BH")
names(ma_p_adj)<-srv_dtMNM_keep$protID

write.csv(ma_p_adj, "busted_output.OUfg_fullReproSet.model_averaged_pvals.csv", quote=TRUE, row.names=TRUE)

#Get proportion of genes under positive selection for tissues and cell types
pos_selec_genes<-names(ma_p_adj)[which(ma_p_adj < 0.05)]
print(paste("Number of genes in analysis:", length(ma_p_adj)))
print(paste("Number of genes with evidence for positive selection:", length(pos_selec_genes)))
prop_selec_all<-length(pos_selec_genes) / length(ma_p_adj)
print(prop_selec_all)

#Read in gene sets for tissues/cell types
som<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

ts<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.testisSpecific.txt", what=character())

bd<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BD",tissue_stringency,"txt", sep="."), what=character())
bu<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BU",tissue_stringency,"txt", sep="."), what=character())
cg<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.CG",tissue_stringency,"txt", sep="."), what=character())
dp<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.DP",tissue_stringency,"txt", sep="."), what=character())
vp<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.VP",tissue_stringency,"txt", sep="."), what=character())
sv<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.SV",tissue_stringency,"txt", sep="."), what=character())

groups<-c("som", "spg", "pre", "spc", "spd", "elo", "bd", "bu", "cg", "dp", "vp", "sv", "ts")
ns<-c(length(ma_p_adj))
prop_selec<-c(prop_selec_all)
prop_test_ps<-c(NA)
for(i in groups){
        print(paste("Working on", i))
        test_set<-ma_p_adj[which(names(ma_p_adj) %in% get(i))]
        print(length(test_set))
        ns<-c(ns, length(test_set))
        print(paste("Number of genes in cell type with evidence for positive selection:", length(which(test_set < 0.05))))
        prop<-length(which(test_set < 0.05)) / length(test_set)
        prop_selec<-c(prop_selec, prop)
        result<-prop.test(length(which(test_set < 0.05)), length(test_set), prop_selec_all)
        prop_test_ps<-c(prop_test_ps, result$p.value)
}

print("Proportion of genes with evidence for positive selection:")
names(prop_selec)<-c("all_genes",groups)
print(prop_selec)

print("Proportion test results, compared to genome-wide positive selection:")
names(prop_test_ps)<-c("all_genes",groups)
print(prop_test_ps)
prop_test_ps_adj<-p.adjust(prop_test_ps, method="BH")
names(prop_test_ps_adj)<-c("all_genes",groups)
print(prop_test_ps_adj)

final_df<-as.data.frame(rbind(ns, prop_selec, prop_test_ps, prop_test_ps_adj))
colnames(final_df)<-c("all_genes",groups)
rownames(final_df)<-c("num_genes", "proportion under pos selec", "pos selec pval", "pos selec BH pval")
write.csv(final_df, paste("selection_results.busted.model_averaged_pval.tissues",tissue_stringency,"csv", sep="."), row.names=TRUE, quote=FALSE)

#Compare testis-specific genes with evidence for selection under each model
srv_dtMNM_padj<-p.adjust(sigdf$srv_dtMNM_pval, method="BH")
names(srv_dtMNM_padj)<-srv_dtMNM_keep$protID
ts_sig_srv_dtMNM_padj<-names(srv_dtMNM_padj)[intersect(which(srv_dtMNM_padj < 0.05), which(names(srv_dtMNM_padj) %in% ts))]

srv_noMNM_padj<-p.adjust(sigdf$srv_noMNM_pval, method="BH")
names(srv_noMNM_padj)<-srv_noMNM_keep$protID
ts_sig_srv_noMNM_padj<-names(srv_noMNM_padj)[intersect(which(srv_noMNM_padj < 0.05), which(names(srv_noMNM_padj) %in% ts))]

noSRV_noMNM_padj<-p.adjust(sigdf$noSRV_noMNM_pval, method="BH")
names(noSRV_noMNM_padj)<-noSRV_noMNM_keep$protID
ts_sig_noSRV_noMNM_padj<-names(noSRV_noMNM_padj)[intersect(which(noSRV_noMNM_padj < 0.05), which(names(noSRV_noMNM_padj) %in% ts))]

ts_sig_map<-intersect(pos_selec_genes, ts)

grps<-c("ts_sig_srv_dtMNM_padj","ts_sig_srv_noMNM_padj","ts_sig_noSRV_noMNM_padj","ts_sig_map")
for(i in 1:(length(grps)-1)){
	for(j in (i+1):length(grps)){
		g1<-grps[i]
		g2<-grps[j]
		print(paste(g1,"vs",g2))
		print(paste("Number of significant testis-specific genes in", g1, length(get(g1))))
		print(paste("Number of significant testis-specific genes in", g2, length(get(g2))))
		print(paste("Number of overlapping significant ts genes between the two models", length(intersect(get(g1), get(g2)))))
	}
}
