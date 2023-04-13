#PURPOSE: Plot proportion of genes under positive selection based on m1a vs m2a in PAML for reproductive tissues
#Also plot propotion of genes that are testis specific 

library(ggplot2)
library(gridExtra)

#Read in and merge dN/dS
posSelec<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.m1aVm2a.txt", header=FALSE)
colnames(posSelec)<-c("protID", "lnL_m2a", "lnL_m1a", "LRT", "raw_p")

#Get subset of genes that converged/ran properly in m1a vs m2a test by removing ones with higher likelihood in m1a model
posSelec_filtered<-c()
for(i in 1:nrow(posSelec)){
	if(!(is.na(posSelec$LRT[i]))){
        	if( (posSelec$lnL_m1a[i]!="NaN") & (posSelec$lnL_m2a[i]!="NaN") ){
                	if(posSelec$lnL_m1a[i] <= posSelec$lnL_m2a[i]){
				posSelec_filtered<-rbind(posSelec_filtered, posSelec[i,])
			}
		}
	}
}

print("Before filtering:")
print(dim(posSelec))
print(head(posSelec))
print("After filtering:")
print(dim(posSelec_filtered))
print(head(posSelec_filtered))

#posSelec_filtered<-posSelec
posSelec_filtered$p.adj<-p.adjust(posSelec_filtered$raw_p, method="BH")
pos_selec_genes<-posSelec_filtered$protID[which(posSelec_filtered$p.adj < 0.05)]
all_selec_genes<-posSelec_filtered$protID
print(paste("Number of genes in analysis:", length(all_selec_genes)))
print(paste("Number of genes with evidence for positive selection:", length(pos_selec_genes)))
prop_selec_all<-length(pos_selec_genes) / length(all_selec_genes)

#Read in gene sets for tissues/cell types
bd<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BD.strict.txt", what=character())
bu<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BU.strict.txt", what=character())
cg<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.CG.strict.txt", what=character())
dp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.DP.strict.txt", what=character())
vp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.VP.strict.txt", what=character())
sv<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.SV.strict.txt", what=character())

#Get proportions and plot
groups<-c("bd", "bu", "cg", "dp", "vp", "sv")
prop_selec<-c()
prop_test_ps<-c()
for(i in groups){
	print(paste("Working on", i))
	test_set<-posSelec_filtered[which(posSelec_filtered$protID %in% get(i)),]
	print(dim(test_set))
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

pdf("prop_pos_selec.m1aVm2a.dean_tissues.pdf")
print(p)
dev.off()

