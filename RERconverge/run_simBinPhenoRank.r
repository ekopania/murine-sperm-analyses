#PURPOSE: Copied code from some of the permulation wrapper scripts to run simBinPhenoRank, instead of simBinPhenoCC

library(RERconverge)
myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/"

source(paste0(myPath,"customPermulationScripts.r"))

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Enter the iteration number for this run")
}
itn<-as.numeric(as.character(args[1]))

myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/"

print("Loading data...")
load(paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds"))

print("Setting foreground species, root, and master tree...")
fg<-c("Hyomys_goliath_ABTC42697", "Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Pseudomys_shortridgei_Z25113", "Notomys_alexis_U1308","Notomys_fuscus_M22830","Notomys_mitchellii_M21518","Zyzomys_pedunculatus_Z34925", "Eropeplus_canus_NMVZ21733", "Paruromys_dominator_JAE4870", "Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074")
sisters_fg<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Eropeplus_canus_NMVZ21733","Paruromys_dominator_JAE4870"), "clade5"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"))
myRoot<-c("Lophuromys_woosnami_LSUMZ37793", "Lophiomys_imhausi_UM5152")
mt<-myTrees$masterTree
#Get trees for all loci
#RERtrees<-myTrees
#Read in trees for a subset of loci
RERtrees<-readTrees(paste0(myPath,"trees_for_RERconverge.noBlanks.tree"), max.read=611, masterTree=mt)
nperms<-100

####RANK PERMULATIONS####
permRank<-getPermsBinaryCustom(numperms=10, fg_vec=fg, sisters_list=sisters_fg, root_sp=myRoot, RERmat=myRER, trees=myTrees, mastertree=mt, permmode="rank")
pdf("perm_trees.pdf", onefile=TRUE, height=15, width=15)
output.binphens<-list()
for(i in 1:nperms){
       perm<-simBinPhenoRank(trees=myTrees, root_sp=myRoot, fg_vec=fg, sisters_list=sisters_fg, plotTreeBool=FALSE)
       plotTreeHighlightBranches(perm, hlspecies=which(perm$edge.length==1), hlcols="red")
       output.binphens[[length(output.binphens)+1]]<-perm
       #print(output.binphens)
}
dev.off()
permulated.binphens<-list()
permulated.binphens[[1]]<-output.binphens
print(length(permulated.binphens[[1]]))
print(permulated.binphens[[1]])
#COPIED FROM PermulationFuncs.R
permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
print(permulated.fg)
#permulated.fg.list = as.list(data.frame(permulated.fg))
phenvec.table = mapply(foreground2Paths,permulated.fg,MoreArgs=list(treesObj=myTrees,clade="all"))
phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])

print("Calculating correlations")
corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=myRER)

print("Generating output permulation table")
#make enrich list/matrices to fill
permPvals=data.frame(matrix(ncol=nperms, nrow=nrow(myRER)))
rownames(permPvals)=rownames(myRER)
permRhovals=data.frame(matrix(ncol=nperms, nrow=nrow(myRER)))
rownames(permRhovals)=rownames(myRER)
permStatvals=data.frame(matrix(ncol=nperms, nrow=nrow(myRER)))
rownames(permStatvals)=rownames(myRER)

for (i in 1:length(corMatList)){
       permPvals[,i] = corMatList[[i]]$P
       permRhovals[,i] = corMatList[[i]]$Rho
       permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
}

permRank=vector("list", 3)
permRank[[1]]=permPvals
permRank[[2]]=permRhovals
permRank[[3]]=permStatvals
names(permRank)=c("corP", "corRho", "corStat")

save(permRank, file=paste0(myPath, paste0("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSimBinPhenoRank.perm", numperms, ".rds")))
