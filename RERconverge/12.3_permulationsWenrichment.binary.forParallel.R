#PURPOSE: Run permulations for binary trait with enrichment

library(RERconverge)
library(biomaRt)
source("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/customPermulationScripts.r")

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Enter the iteration number for this run")
}
itn<-as.numeric(as.character(args[1]))

myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/"

print("Loading data...")
load(paste0(myPath,"RERconverge_output.logRTM_binary.enrichment.rds"))
load(paste0(myPath,"RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds"))

print("Setting foreground species, root, and master tree...")
fg<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074")
#fg_sisters<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308", "Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074"))
myRoot<-c("Lophuromys_woosnami_LSUMZ37793", "Lophiomys_imhausi_UM5152")
#myRoot<-NULL
RERtrees<-myTrees
mt<-RERtrees$masterTree

print("Loading annotations and changing mouse protein IDs to human gene symbols...")
#load annotations
annots<-read.gmt("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/c2.all.v7.4.symbols.gmt")
annotlist=list(annots)
names(annotlist)="MSigDBpathways"

##Change mouse protein IDs to human gene symbols
###Doing res and RER matrix at the same time based on the same protein IDs, so make sure they're the same and in the same order
#stopifnot(all.equal(rownames(myRER), rownames(res)))
#mouse_prot<-unlist(sapply(rownames(myRER), function(x) gsub("-.*","",x)))
###Using an old version of Ensembl because getLDS() has issues after the April 2022 update
#ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org/")
#ens_hum<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "https://nov2020.archive.ensembl.org/")
#all_genes<-getBM(attributes=c('ensembl_peptide_id','external_gene_name'), mart=ens_mus)
#colnames(all_genes)<-c("prot_id","gene_name")
#musNames<-unlist(sapply(mouse_prot, function(x) all_genes$gene_name[which(all_genes$prot_id==x)]))
#hum_genes<-getLDS(attributes=c("ensembl_peptide_id","external_gene_name"), values=musNames, mart=ens_mus, attributesL=c("external_gene_name"), martL=ens_hum, uniqueRows=T)
#colnames(hum_genes)<-c("mus_ID", "mus_name", "hum_name")
#print(dim(hum_genes))
#print(head(hum_genes))
#rownames(myRER)<-mouse_prot
#print(myRER[1:5, 1:5])
#rownames(res)<-mouse_prot
#print(head(res))
#
###Some stuff to deal with 1:many and many:1 orthos, as well as genes that don't have an ortholog in both species
###If many human genes, keep the first one
###If many mouse genes, keep the first one
###If no ortholog, remove from dataset
#hum_genes_sorted<-hum_genes[order(hum_genes$mus_ID),]
#hum_genes_uniq<-hum_genes_sorted[which(!(duplicated(hum_genes_sorted$mus_ID))),]
#hum_genes_final<-hum_genes_uniq[which(hum_genes_uniq$mus_ID != ""),]
#
#myRER_humNames<-c()
#res_humNames<-c()
#no_ortho<-c()
#multi_mouse<-c()
#for(i in 1:nrow(myRER)){
#       if(rownames(myRER)[i] %in% hum_genes_final$mus_ID){
#               hum_g<-hum_genes_final$hum_name[which(hum_genes_final$mus_ID==rownames(myRER)[i])]
#                       if(hum_g %in% rownames(myRER_humNames)){
#                               #print(paste("Two mouse genes for the following human gene:", hum_g))
#                               multi_mouse<-c(multi_mouse, rownames(myRER)[i])
#                               #Don't add because then there will be multiple rows with the same name (same human gene) in the RERmatrix, which will give an error
#                               #myRER_humNames<-rbind(myRER_humNames, myRER[i,])
#                               #rownames(myRER_humNames)[nrow(myRER_humNames)]<-hum_g
#                       } else{
#                               myRER_humNames<-rbind(myRER_humNames, myRER[i,])
#                               rownames(myRER_humNames)[nrow(myRER_humNames)]<-hum_g
#                               res_humNames<-rbind(res_humNames, res[i,])
#                               rownames(res_humNames)[nrow(res_humNames)]<-hum_g
#                       }
#       } else{
#               no_ortho<-c(no_ortho, rownames(myRER)[i])
#       }
#}
#print(paste("Number of genes removed from myRER due to no ortholog:", length(no_ortho)))
#print(paste("Number of genes with multiple mouse genes corresponding to a single human ortholog:", length(multi_mouse)))
#print(dim(myRER_humNames))
#print(myRER_humNames[1:5,1:5])
#save(myRER_humNames, res_humNames, file="RERconverge_output.logRTM_binary.permWenrichment.rds")

load(paste0(myPath,"submit/RERconverge_output.logRTM_binary.permWenrichment.rds"))

#perform permulations - complete case with enrichment
print("Running permulations...")
#permCCwEnrich<-getPermsBinary(numperms=10, fg_vec=fg, sisters_list=fg_sisters, root_sp=myRoot, RERmat=myRER_humNames, trees=myTrees, mastertree=mt, permmode="cc",calculateenrich=TRUE, annotlist=annotlist, min.pos=4)
#RERtrees<-RERtrees$trees[1:10] #First 10 trees for testing
permSSMwEnrich<-getPermsBinaryCustom(numperms=2, fg_vec=fg, sisters_list=NULL, root_sp=myRoot, RERmat=myRER_humNames, trees=RERtrees, mastertree=mt, permmode="ssm", calculateenrich=TRUE, annotlist=annotlist)
#save(myRER_humNames, res_humNames, permCCwEnrich, file="RERconverge_output.logRTM_binary.permWenrichment.rds")
save(myRER_humNames, res_humNames, permSSMwEnrich, file=paste0(myPath,"RERconverge_output.logRTM_binary.permWenrichment.parallel",itn,".rds"))
q()

print("Permulation results for enrichment:")
#Get p-value and correct for multiple tests
#enrichpermpvals<-permpvalenrich(enrichment, permCCwEnrich)
enrichpermpvals<-permpvalenrich(enrichment, permSSMwEnrich)
enrichpermpvals_adj<-p.adjust(enrichpermpvals[[1]], method="BH")

repro_sets<-c("MATZUK_SPERMATID_DIFFERENTIATION", "MATZUK_SPERMATOCYTE", "MATZUK_SPERMATOGONIA", "MATZUK_SPERMATOZOA", "REACTOME_SPERM_MOTILITY_AND_TAXES", "WEBER_METHYLATED_HCP_IN_SPERM_DN", "WEBER_METHYLATED_HCP_IN_SPERM_UP", "WEBER_METHYLATED_ICP_IN_SPERM_DN", "WEBER_METHYLATED_ICP_IN_SPERM_UP", "WEBER_METHYLATED_LCP_IN_SPERM_DN", "WEBER_METHYLATED_LCP_IN_SPERM_UP", "CHEN_ETV5_TARGETS_TESTIS", "CHEN_ETV5_TARGETS_SERTOLI", "REACTOME_TRANSCRIPTIONAL_REGULATION
+ _OF_TESTIS_DIFFERENTIATION", "SU_TESTIS", "YOKOE_CANCER_TESTIS_ANTIGENS", "MATZUK_MALE_REPRODUCTION_SERTOLI", "MATZUK_CENTRAL_FOR_FEMALE_FERTILITY", "MATZUK_CUMULUS_EXPANSION", "MATZUK_EARLY_ANTRAL_FOLLICLE", "MATZUK_FERTILIZATION", "MATZUK_IMPLANTATION_AND_UTERINE", "MATZUK_LUTEAL_GENES", "MATZUK_MEIOTIC_AND_DNA_REPAIR", "MATZUK_OVULATION", "MATZUK_POSTIMPLANTATION_AND_POSTPAR
+ TUM", "MATZUK_PREOVULATORY_FOLLICLE","MATZUK_STEROIDOGENESIS","REACTOME_REPRODUCTION", "WP_MALE_INFERTILITY")

print(paste("Number of repro gene sets in permulated enrichment result:", length(which(repro_sets %in% names(enrichpermpvals_adj)))))
print("Enrichment results for repro sets:")
print(enrichpermpvals_adj[repro_sets])

print("Getting permulation results for res; saving to a table...")
#corpermpvals<-permpvalcor(res_humNames, permCCwEnrich)
corpermpvals<-permpvalcor(res_humNames, permSSMwEnrich)
res_humNames$permpval<-corpermpvals[match(rownames(res), names(corpermpvals))]
res_humNames$permpvaladj<-p.adjust(res$permpval, method="BH")
count=1
while(count<=length(enrichment)){
	enrichment[[count]]$permpval=enrichpermpvals[[count]][match(rownames(enrichment[[count]]),
	names(enrichpermpvals[[count]]))]
	enrichment[[count]]$permpvaladj=p.adjust(enrichment[[count]]$permpval, method="BH")
	count=count+1
}
enrichment_sorted<-enrichment$MSigDBpathways[order(enrichment$MSigDBpathways$permpvaladj),]

#save data
#save(myRER_humNames, res_humNames, permCCwEnrich, enrichpermpvals, enrichpermpvals_adj, enrichment, enrichment_sorted, file="RERconverge_output.logRTM_binary.permWenrichment.rds")
save(myRER_humNames, res_humNames, permSSMwEnrich, enrichpermpvals, enrichpermpvals_adj, enrichment, enrichment_sorted, file="RERconverge_output.logRTM_binary.permWenrichment.rds")
write.table(res, file="RERconverge_output.logRTM_binary_residuals.permWenrichment.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(enrichment_sorted, file="RERconverge_output.logRTM_binary.enrichments.permulated.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)

print("Done with 12.2_permulationsWenrichment.binary.R")
