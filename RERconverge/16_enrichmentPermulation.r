#PURPOSE: Run enrichment permulations

library(topGO)
library(RERconverge)
source("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/GO_permulation.r")

#Load data
myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/"

print("Loading data...")
load(paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds"))

fg<-c("Hyomys_goliath_ABTC42697", "Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Pseudomys_shortridgei_Z25113", "Notomys_alexis_U1308","Notomys_fuscus_M22830","Notomys_mitchellii_M21518","Zyzomys_pedunculatus_Z34925", "Eropeplus_canus_NMVZ21733", "Paruromys_dominator_JAE4870", "Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074")
sisters_fg<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Eropeplus_canus_NMVZ21733","Paruromys_dominator_JAE4870"), "clade5"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"))
#fg<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308","Notomys_fuscus_M22830","Notomys_mitchellii_M21518","Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074")
#sisters_fg<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"))
print("Foregrounds:")
print(fg)

myRoot<-c("Lophuromys_woosnami_LSUMZ37793", "Lophiomys_imhausi_UM5152")
nperms<-1000

load(paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.customSisters.rds"))

#Format RER matrix
prots<-unlist(sapply(rownames(myRER), function(x) gsub("-.*","",x)))
genes<-protToGene(prots)
print(dim(myRER))
print(length(genes))
print(head(genes))
rownames(myRER)<-genes

#Format permulation matrix
stopifnot(all.equal(rownames(combpermSSM$corP), rownames(combpermSSM$corRho)))
stopifnot(all.equal(rownames(combpermSSM$corP), rownames(combpermSSM$corStat)))
prots<-unlist(sapply(rownames(combpermSSM$corP), function(x) gsub("-.*","",x)))
genes<-protToGene(prots)
dim(combpermSSM$corP)
print(length(genes))
print(head(genes))
na.count<-1
#Number the NAs so they are unique and can be appended as rownames
for(g in 1:length(genes)){
	if(is.na(genes[g])){
		genes[g]<-paste("NA",na.count,sep=".")
		na.count<-na.count+1
	}
}
rownames(combpermSSM$corP)<-genes
rownames(combpermSSM$corRho)<-genes
rownames(combpermSSM$corStat)<-genes

#Run enrichment permulations
print("Running enrichment perms...")
enrich_perms<-getGOEnrichPermsBinary(numperms=nperms, fg_vec=fg, sisters_list=sisters_fg, root_sp=myRoot, RERmat=myRER, trees=myTrees, method="k", min.pos=2, combpermSSM, GO_thresh=250)

#Save output
print("Saving output...")
save(enrich_perms, file=paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.enrichment.customSisters.rds"))
