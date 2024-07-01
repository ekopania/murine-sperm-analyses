#PURPOSE: Test for GO enrichment in top n% of genes from RERconverge output

library(RERconverge)
library(topGO)
library(biomaRt)

geneset<-"top500" #top100, top500, top1000, p05, p01, p11
dataset<-"_binary.q90.RTMspeciesOnly" #"_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.fromMasterTree" #"", "_binary.RTMspeciesOnly", "_binary.q10.RTMspeciesOnly", "_binary.q20.RTMspeciesOnly", "_binary.q80.RTMspeciesOnly", "_binary.q90.RTMspeciesOnly", "_binary.q50low.RTMspeciesOnly", "_binary.q50high.RTMspeciesOnly", ""_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.customSisters"
print(paste("Running GO enrichment analysis for gene set", geneset, "and dataset", dataset))

#Read in data
print("Loading data...")
load(paste0("RERconverge_output.logRTM", dataset, ".rds"))
#With permulations:
#res<-read.csv("RERconverge_output.logRTM_binary_residuals.OUmodel.RTMspeciesOnly.withPermulationsSSM.customSSM.1000perms.csv", row.names=1)
#res<-read.csv("RERconverge_output.logRTM_binary_residuals.OUmodel.RTMspeciesOnly.withPermulationsSSM.customSSM.fromMasterTree.csv", row.names=1)

print(head(res))

#Remove paralogs from data
paralogs<-scan("paralog_loci.txt", what=character())
res<-res[which(!(rownames(res) %in% paralogs)),]

#Generate histograms of various RERconverge results parameters; use these to determiine cutoffs for GO enrichment
#print("Generating histograms of RERconverge results...")
#pdf("RERconverge_output.logRTM_binary.results_histograms.pdf", onefile=TRUE)
#hist(res$Rho, main="Rho value histogram")
#hist(res$P, main="P value histogram")
#hist(res$p.adj, main="Adjusted P value histogram")
#signP<-c()
#for(i in 1:nrow(res)){
#	if(is.na(res$Rho[i])){
#		signP<-c(signP, NA)
#	} else{
#		if(res$Rho[i] < 0){
#			signP<-c(signP, 0-res$P[i])
#		} else{
#			signP<-c(signP, res$P[i])
#		}
#	}
#}
#hist(signP, main="Signed P-value based on sign of Rho")
#dev.off()
#Looks like uncorrected p < 0.05 will be a good starting point (2022-04-20)

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
#print(length(unique(sort(all_genes_inRes$prot_id))))
#print(length(unique(sort(all_genes_inRes$gene_id))))
#print(length(unique(sort(all_genes_inRes$gene_name))))
res_hasGeneID<-res[which(rownames(res) %in% all_genes_inRes$prot_id),]
stopifnot(all.equal(rownames(res_hasGeneID), all_genes_inRes$prot_id))
mouse_gene<-unlist(sapply(rownames(res_hasGeneID), function(x) all_genes_inRes$gene_id[which(all_genes_inRes$prot_id==x)]))
stopifnot(length(mouse_gene)==nrow(res_hasGeneID))
rownames(res_hasGeneID)<-mouse_gene
print(paste("Number of genes lost in converting to gene ID:", nrow(res)-nrow(res_hasGeneID)))

print("Running topGO analysis...")
#CHANGE "subset" for different threshold of top or significant genes
if(geneset=="p05"){
	#Uncorrected P-value < 0.05
	subset<-rownames(res_hasGeneID)[which(res_hasGeneID$permpval < 0.05)]
} else if(geneset=="p01"){
	#Uncorrected P-value < 0.01
	subset<-rownames(res_hasGeneID)[which(res_hasGeneID$permpval < 0.01)]
} else if(geneset=="p11"){
        #Corrected P-value < 0.11
        #subset<-rownames(res_hasGeneID)[which(res_hasGeneID$P < 0.01)]
	subset<-rownames(res_hasGeneID)[which(res_hasGeneID$permpval < 0.11)]
} else if(geneset=="top1000"){
        #Top 500 genes
        #subset<-rownames(res_hasGeneID[order(res_hasGeneID$P),])[1:500]
        subset<-rownames(res_hasGeneID[order(res_hasGeneID$permpval),])[1:1000]
} else if(geneset=="top500"){
	#Top 500 genes
	subset<-rownames(res_hasGeneID[order(res_hasGeneID$P),])[1:500]
	#subset<-rownames(res_hasGeneID[order(res_hasGeneID$permpval),])[1:500]
} else if(geneset=="top100"){
	#Top 100 genes
        #subset<-rownames(res_hasGeneID[order(res_hasGeneID$P),])[1:100]
	subset<-rownames(res_hasGeneID[order(res_hasGeneID$permpval),])[1:100]
} else{
	stop("Invalid gene set")
}
print(paste("Number of genes in focal subset:", length(subset)))
print(paste("Total number of genes in dataset:", nrow(res_hasGeneID)))
total<-rownames(res_hasGeneID)
inSubset<-c()
for(gene in total){
	if(gene %in% subset){
		inSubset<-c(inSubset, 1)
	} else{
		inSubset<-c(inSubset, 0)
	}
}
if(length(subset) != length(which(inSubset==1))){
	print("ERROR: length of 'subset' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
	return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("RERconverge_output.logRTM",dataset,".GOenrichment.",geneset,".csv"), row.names=FALSE)

#Repeat for accelerated and decerlated separately
res_acc<-res_hasGeneID[which(res_hasGeneID$Rho > 0), ]
res_dec<-res_hasGeneID[which(res_hasGeneID$Rho <= 0), ]
#res_acc<-res_hasGeneID[which(res_hasGeneID$permstats > 0), ]
#res_dec<-res_hasGeneID[which(res_hasGeneID$permstats <= 0), ]
if(geneset=="p05"){
        #Uncorrected P-value < 0.05
        subset_acc<-rownames(res_acc)[which(res_acc$permpval < 0.05)]
	subset_dec<-rownames(res_dec)[which(res_dec$permpval < 0.05)]
} else if(geneset=="p01"){
        #Uncorrected P-value < 0.01
        subset_acc<-rownames(res_acc)[which(res_acc$permpval < 0.01)]
	subset_dec<-rownames(res_dec)[which(res_dec$permpval < 0.01)]
} else if(geneset=="p11"){
        #Corrected P-value < 0.11
        #subset_acc<-rownames(res_acc)[which(res_acc$P < 0.01)]
        #subset_dec<-rownames(res_dec)[which(res_dec$P < 0.01)]
	subset_acc<-rownames(res_acc)[which(res_acc$permpval < 0.11)]
        subset_dec<-rownames(res_dec)[which(res_dec$permpval < 0.11)]
} else if(geneset=="top1000"){
        #Top 500 genes
        #subset_acc<-rownames(res_acc[order(res_acc$P),])[1:250]
        #subset_dec<-rownames(res_dec[order(res_dec$P),])[1:250]
        subset_acc<-rownames(res_acc[order(res_acc$permpval),])[1:500]
        subset_dec<-rownames(res_dec[order(res_dec$permpval),])[1:500]
} else if(geneset=="top500"){
        #Top 500 genes
        subset_acc<-rownames(res_acc[order(res_acc$P),])[1:250]
	subset_dec<-rownames(res_dec[order(res_dec$P),])[1:250]
	#subset_acc<-rownames(res_acc[order(res_acc$permpval),])[1:250]
        #subset_dec<-rownames(res_dec[order(res_dec$permpval),])[1:250]
} else if(geneset=="top100"){
        #Top 100 genes
        #subset_acc<-rownames(res_acc[order(res_acc$P),])[1:50]
	#subset_dec<-rownames(res_dec[order(res_dec$P),])[1:50]
	subset_acc<-rownames(res_acc[order(res_acc$permpval),])[1:50]
        subset_dec<-rownames(res_dec[order(res_dec$permpval),])[1:50]
} else{
        stop("Invalid gene set")
}

print(paste("Number of genes in focal accelerated subset:", length(subset_acc)))
print(paste("Total number of genes in accelerated dataset:", nrow(res_acc)))
print(head(subset_acc))
total<-rownames(res_acc)
inSubset<-c()
for(gene in total){
        if(gene %in% subset_acc){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_acc) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_acc' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_acc", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("RERconverge_output.logRTM",dataset,".GOenrichment.acceleratedONLY.",geneset,".csv"), row.names=FALSE)

print(paste("Number of genes in focal decelerated subset:", length(subset_dec)))
print(paste("Total number of genes in decelerated dataset:", nrow(res_dec)))
total<-rownames(res_dec)
inSubset<-c()
for(gene in total){
        if(gene %in% subset_dec){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_dec) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_dec' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_dec", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("RERconverge_output.logRTM",dataset,".GOenrichment.deceleratedONLY.",geneset,".csv"), row.names=FALSE)

print("Done with 13_GOenrichment.r")
