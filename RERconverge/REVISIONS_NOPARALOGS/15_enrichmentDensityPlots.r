#PURPOSE: Generate plots of enrichment density by ranked P-value using scripts from Amanda and Maria

library(topGO)
library(biomaRt)
source("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/MariaBarcodePlotFinal.R")

paralogs<-scan("paralog_loci.txt", what=character())

print("Reading in and formatting data...")
rhoMat<-c()
signedPMat<-c()
#Loop through all binary RERconverge runs for RTM
#RERrun<-c("_binary.OUmodel.RTMspeciesOnly", "_binary.q10.RTMspeciesOnly", "_binary.q20.RTMspeciesOnly", "_binary.q80.RTMspeciesOnly", "_binary.q90.RTMspeciesOnly", "")
#Just lower RTM and continuous
#RERrun<-c("_binary.OUmodel.RTMspeciesOnly", "_binary.q10.RTMspeciesOnly", "_binary.q20.RTMspeciesOnly", "_binary.q50low.RTMspeciesOnly", "")
#Just higher RTM and continuous
#RERrun<-c("_binary.q80.RTMspeciesOnly", "_binary.q90.RTMspeciesOnly", "_binary.q50high.RTMspeciesOnly", "")
#Just OU model
#RERrun<-c("_binary.OUmodel.RTMspeciesOnly")
#OU model, continuous, and top 10% RTM
RERrun<-c("_binary.OUmodel.RTMspeciesOnly", "_binary.q90.RTMspeciesOnly", "_continuous.RTMspeciesOnly")
for(i in RERrun){
	load(paste0("RERconverge_output.logRTM", i, ".rds"))
	res<-res[which(!(rownames(res) %in% paralogs)),]
	print(head(res))
	print("Converting mouse protein IDs to gene names...")
	mouse_prot<-unlist(sapply(rownames(res), function(x) gsub("-.*","",x)))
	rownames(res)<-mouse_prot
	res<-res[order(rownames(res)),]
	ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org/")
	all_genes<-getBM(attributes=c('ensembl_peptide_id','ensembl_gene_id','external_gene_name'), mart=ens_mus)
	colnames(all_genes)<-c("prot_id","gene_id","gene_name")
	all_genes_sorted<-all_genes[order(all_genes$prot_id),]
	all_genes_uniq<-all_genes_sorted[which(!(duplicated(all_genes_sorted$prot_id))),]
	all_genes_final<-all_genes_uniq[which(all_genes_uniq$prot_id != ""),]
	all_genes_inRes<-all_genes_final[which(all_genes_final$prot_id %in% rownames(res)),]
	print(dim(all_genes_inRes))
	print(head(all_genes_inRes))
	print(tail(all_genes_inRes))
	res_hasGeneID<-res[which(rownames(res) %in% all_genes_inRes$prot_id),]
	print(dim(res_hasGeneID))
        print(head(res_hasGeneID))
        print(tail(res_hasGeneID))
	stopifnot(all.equal(rownames(res_hasGeneID), all_genes_inRes$prot_id))
	mouse_gene<-unlist(sapply(rownames(res_hasGeneID), function(x) all_genes_inRes$gene_id[which(all_genes_inRes$prot_id==x)]))
	stopifnot(length(mouse_gene)==nrow(res_hasGeneID))
	rownames(res_hasGeneID)<-mouse_gene
	print(paste("Number of genes lost in converting to gene ID:", nrow(res)-nrow(res_hasGeneID)))
	#Make sure gene names the same and in the same order across different datasets
	if(is.vector(rhoMat)){
		stopifnot(all.equal(names(rhoMat), rownames(res_hasGeneID)))
	} else if(length(rhoMat) > 0){
		stopifnot(all.equal(rownames(rhoMat), rownames(res_hasGeneID)))
	}
	#Make the matrix such that pos rho always associated w/ lower RTM; add gene names as rownames
	if(i %in% c("_binary.OUmodel.RTMspeciesOnly", "_binary.q10.RTMspeciesOnly", "_binary.q20.RTMspeciesOnly", "_binary.q50low.RTMspeciesOnly")){
		rhoMat<-cbind(rhoMat, res_hasGeneID$Rho)
	} else{
		rhoMat<-cbind(rhoMat, -(res_hasGeneID$Rho))
	}
	if(is.vector(rhoMat)){
		names(rhoMat)<-rownames(res_hasGeneID)
	} else{
		rownames(rhoMat)<-rownames(res_hasGeneID)
	}
}

print("rhoMat info:")
#colnames(rhoMat)<-c("OUmodel","low10","low20","high20","high10", "continuous")
#colnames(rhoMat)<-c("OUmodel","low10","low20","low50","continuous")
#colnames(rhoMat)<-c("high20","high10", "high50", "continuous")
#colnames(rhoMat)<-c("OUmodel")
colnames(rhoMat)<-c("OUmodel","high10", "continuous")
print(paste("Matrix?", is.matrix(rhoMat)))
print(dim(rhoMat))
print(head(rhoMat))

#For signed P, will have to use apply to get the sign of rho and apply it to p value; will deal with that later

print("Getting gene sets...")
annots<-annFUN.org("BP", mapping="org.Mm.eg.db", ID="ensembl")
GO_spermatogenesis<-unlist(annots["GO:0007283"])
#print(head(GO_spermatogenesis))
GO_maleGameteGeneration<-unlist(annots["GO:0048232"])
#print(head(GO_maleGameteGeneration))
GO_gameteGeneration<-unlist(annots["GO:0007276"])
#Can add others
##Read in lists of genes enriched in different spermatogenesis cell types or male repro tissues
LZ_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt", what=character())
RS_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt", what=character())
TS_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_testisSpecific.txt", what=character())
som_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spt_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())
bd_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BD.any.txt", what=character())
bu_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BU.any.txt", what=character())
cg_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.CG.any.txt", what=character())
dp_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.DP.any.txt", what=character())
vp_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.VP.any.txt", what=character())
sv_prot<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.SV.any.txt", what=character())
LZ_genes<-unique(sort(unlist(sapply(LZ_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
RS_genes<-unique(sort(unlist(sapply(RS_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
TS_genes<-unique(sort(unlist(sapply(TS_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
som_genes<-unique(sort(unlist(sapply(som_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
spg_genes<-unique(sort(unlist(sapply(spg_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
pre_genes<-unique(sort(unlist(sapply(pre_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
spc_genes<-unique(sort(unlist(sapply(spc_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
spt_genes<-unique(sort(unlist(sapply(spt_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
elo_genes<-unique(sort(unlist(sapply(elo_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
bd_genes<-unique(sort(unlist(sapply(bd_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
bu_genes<-unique(sort(unlist(sapply(bu_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
cg_genes<-unique(sort(unlist(sapply(cg_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
dp_genes<-unique(sort(unlist(sapply(dp_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
vp_genes<-unique(sort(unlist(sapply(vp_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
sv_genes<-unique(sort(unlist(sapply(sv_prot, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))


print("Running plot function...")
#pdf("enrichment_density.pdf", height=5, width=7, onefile=TRUE)
#pdf("enrichment_density.lowerFGs.span0_1.pdf", height=5, width=7, onefile=TRUE)
#pdf("enrichment_density.higherFGs.span0_1.pdf", height=5, width=7, onefile=TRUE)
#pdf("enrichment_density.OUmodel.span0_2.pdf", height=5, width=7, onefile=TRUE)
pdf("enrichment_density.OUmodel_high10_cont.span0_2.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(rhoMat, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis")
makeMultiGSEAplot(rhoMat, GO_maleGameteGeneration, span=.2, title="GO category: male gamete generation")
makeMultiGSEAplot(rhoMat, GO_gameteGeneration, span=.2, title="GO category: gamete generation")
makeMultiGSEAplot(rhoMat, LZ_genes, span=.2, title="Induced early")
makeMultiGSEAplot(rhoMat, RS_genes, span=.2, title="Induced late")
makeMultiGSEAplot(rhoMat, TS_genes, span=.2, title="Testis-specific")
makeMultiGSEAplot(rhoMat, som_genes, span=.2, title="SC - Somatic")
makeMultiGSEAplot(rhoMat, spg_genes, span=.2, title="SC - Spermatogonia")
makeMultiGSEAplot(rhoMat, pre_genes, span=.2, title="SC - Pre-leptotene")
makeMultiGSEAplot(rhoMat, spc_genes, span=.2, title="SC - Spermatocytes")
makeMultiGSEAplot(rhoMat, spt_genes, span=.2, title="SC - Spermatids")
makeMultiGSEAplot(rhoMat, elo_genes, span=.2, title="SC - Elongating")
makeMultiGSEAplot(rhoMat, bd_genes, span=.2, title="AG - Bulbourethral diverticulum")
makeMultiGSEAplot(rhoMat, bu_genes, span=.2, title="AG - Bulbourethral gland")
makeMultiGSEAplot(rhoMat, cg_genes, span=.2, title="AG - Coagulating gland")
makeMultiGSEAplot(rhoMat, dp_genes, span=.2, title="AG - Dorsal prostate")
makeMultiGSEAplot(rhoMat, vp_genes, span=.2, title="AG - Ventral prostate")
makeMultiGSEAplot(rhoMat, sv_genes, span=.2, title="AG - Seminal vesicle")
dev.off()

q()


print("Plotting with selection results...")
source("MariaBarcodePlotFinal.colorBySelec.R")
relax<-read.table("/mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.sigOnly.posOnly.m1aVbsnull.txt", header=TRUE)
ps<-read.table("/mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.sigOnly.posOnly.m1aVm2a.txt", header=TRUE)
bs<-read.table("/mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.sigOnly.posOnly.txt", header=TRUE)
relax_genes<-unique(sort(unlist(sapply(relax$protID, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
ps_genes<-unique(sort(unlist(sapply(ps$protID, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
bs_genes<-unique(sort(unlist(sapply(bs$protID, function(x) all_genes$gene_id[which(all_genes$prot_id == x)]))))
selec_status<-c()
#1 = relaxed purifying selection; 2 = pervasive positive selection; 3 = branch-site positive selection; 0 = none
for(i in rownames(rhoMat)){
	if(i %in% relax_genes){
		selec_status<-c(selec_status, 1)
	} else if(i %in% ps_genes){
		selec_status<-c(selec_status, 2)
	} else if(i %in% bs_genes){
		selec_status<-c(selec_status, 3)
	} else{
		selec_status<-c(selec_status, 0)
	}
}
rhoMat_selec<-as.matrix(cbind(rhoMat, selec_status))
#rhoMat_selec[,"OUmodel"]<-as.numeric(as.character(rhoMat_selec[,"OUmodel"]))
print(head(rhoMat_selec))

rhoMat_selec_spermOnly<-rhoMat_selec[which(rownames(rhoMat_selec) %in% GO_spermatogenesis),]
print("rhoMat for spermatogenesis genes:")
print(dim(rhoMat_selec_spermOnly))
print(head(rhoMat_selec_spermOnly))
print(paste("Number of genes under relaxed selection:", length(which(rhoMat_selec_spermOnly[,'selec_status']==1))))
print(paste("Number of genes under positive selection:", length(which(rhoMat_selec_spermOnly[,'selec_status']==2))))
print(paste("Number of genes under branch-site positive selection:", length(which(rhoMat_selec_spermOnly[,'selec_status']==3))))
print(paste("Number of genes in none of the selection categories:", length(which(rhoMat_selec_spermOnly[,'selec_status']==0))))

pdf("enrichment_density.OUmodel.span0_2.withSelection.pdf", height=5, width=10, onefile=TRUE)
myCol<-c("#ffffff", "#028090", "#3E4768", "#984447", "#000000")
#none = #CCCCCC, relaxed = #028090, pos selec = #2A2E45, bs pos selec = #984447, worm = #000000
makeMultiGSEAplotSelection(rhoMat_selec, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis", cols=myCol)

#THIS DIDN'T WORK - trying something new...see below
#Color worms by selection status; remove genes in "none" category
#Make a new matrix where each selection category is its own column
#print("Coloring worms by selection...")
#relax_col<-c()
#pos_col<-c()
#bs_col<-c()
#none_col<-c()
#for(i in 1:nrow(rhoMat_selec)){
#	if(rhoMat_selec[i, "selec_status"] == 1){
#		relax_col<-c(relax_col, rhoMat_selec[i,'OUmodel'])
#	} else{
#		relax_col<-c(relax_col, NA)
#	}
#	if(rhoMat_selec[i, "selec_status"] == 2){
#		pos_col<-c(pos_col, rhoMat_selec[i,'OUmodel'])
#        } else{
#                pos_col<-c(pos_col, NA)
#        }
#	if(rhoMat_selec[i, "selec_status"] == 3){
#                bs_col<-c(bs_col, rhoMat_selec[i,'OUmodel'])
#        } else{
#                bs_col<-c(bs_col, NA)
#        }
#	if(rhoMat_selec[i, "selec_status"] == 0){
#                none_col<-c(none_col, rhoMat_selec[i,'OUmodel'])
#        } else{
#                none_col<-c(none_col, NA)
#        }
#}
#rhoMat_bySelec<-as.matrix(cbind(relaxed=relax_col, positive=pos_col, branch_site=bs_col, none=none_col, selec_status=rhoMat_selec[,'selec_status']))
##rhoMat_bySelec<-rhoMat_bySelec_all[which(rhoMat_bySelec_all[,'selec_status'] > 0),]
#print(head(rhoMat_bySelec))
#print(dim(rhoMat_bySelec))
#sperm_selec_genes<-rownames(rhoMat_bySelec)[intersect(which(rhoMat_bySelec[,'selec_status']!=0), which(rownames(rhoMat_bySelec) %in% GO_spermatogenesis))]
#print(length(sperm_selec_genes))
#print(head(sperm_selec_genes))
#myCol<-c("#028090", "#2A2E45", "#984447", "#000000", "#000000", "#984447", "#2A2E45", "#028090") #Same colors as before; goes 1,2,3 and then alphabetically
#makeMultiGSEAplotSelectionWithWorms(rhoMat_bySelec, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; by selection", cols=myCol)
#dev.off()

print("Using genes with both selection information and in spermatogenesis GO category as input genes")
GO_spermatogenesis_relaxed<-intersect(rownames(rhoMat_selec)[which(rhoMat_selec[,'selec_status']==1)], GO_spermatogenesis)
print(length(GO_spermatogenesis_relaxed))
makeMultiGSEAplotSelection(rhoMat_selec, GO_spermatogenesis_relaxed, span=.2, title="GO category: spermatogenesis - relaxed selection", cols=c("#028090","#028090"))
GO_spermatogenesis_positive<-intersect(rownames(rhoMat_selec)[which(rhoMat_selec[,'selec_status']==2)], GO_spermatogenesis)
print(length(GO_spermatogenesis_positive))
makeMultiGSEAplotSelection(rhoMat_selec, GO_spermatogenesis_positive, span=.2, title="GO category: spermatogenesis - positive selection", cols=c("#3E4768","#3E4768"))
GO_spermatogenesis_branch_site<-intersect(rownames(rhoMat_selec)[which(rhoMat_selec[,'selec_status']==3)], GO_spermatogenesis)
print(length(GO_spermatogenesis_branch_site))
makeMultiGSEAplotSelection(rhoMat_selec, GO_spermatogenesis_branch_site, span=.2, title="GO category: spermatogenesis - branch_site selection", cols=c("#984447","#984447"))
dev.off()

rhoMat_selec_ordered<-rhoMat_selec[order(-rhoMat_selec[,'OUmodel']),]
write.table(rhoMat_selec_ordered, file="rhoMatrix_withSelectionCategory.ordered.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
write.table(rhoMat_selec, file="rhoMatrix_withSelectionCategory.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

print("Done with 15_enrichmentDensityPlots.r")
