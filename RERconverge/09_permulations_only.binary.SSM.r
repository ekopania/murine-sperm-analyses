#PURPOSE: Load data from a previous BINARY RERconverge run and run permulations

library(RERconverge)
myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/"

source(paste0(myPath,"customPermulationScripts.r"))

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Enter the iteration number for this run")
}
itn<-as.numeric(as.character(args[1]))

print("Loading data...")
load(paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds"))

print("Setting foreground species, root, and master tree...")
fg<-c("Hyomys_goliath_ABTC42697", "Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Pseudomys_shortridgei_Z25113", "Notomys_alexis_U1308","Notomys_fuscus_M22830","Notomys_mitchellii_M21518","Zyzomys_pedunculatus_Z34925", "Eropeplus_canus_NMVZ21733", "Paruromys_dominator_JAE4870", "Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074")
sisters_fg<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Eropeplus_canus_NMVZ21733","Paruromys_dominator_JAE4870"), "clade5"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"))
myRoot<-c("Lophuromys_woosnami_LSUMZ37793", "Lophiomys_imhausi_UM5152")
mt<-myTrees$masterTree
#Get trees for all loci
RERtrees<-myTrees
#Read in trees for a subset of loci
#RERtrees<-readTrees(paste0(myPath,"trees_for_RERconverge.noBlanks.tree"), max.read=611, masterTree=mt)
nperms<-50

#print("PATHVEC:")
#print(pathvec)

print("Running permulations...")
#Run permulations
#perform binary species subset match (SSM) permulation - WITH PLOTTING
#pdf(paste0("perm_trees.SSMperms",nperms,".parallel",itn,".pdf"), onefile=TRUE, height=15, width=15)
#permSSM<-getPermsBinaryQuiet(numperms=nperms, fg_vec=fg, sisters_list=sisters_fg, root_sp=myRoot, RERmat=myRER, trees=RERtrees, mastertree=mt, permmode="ssm", trees_list=myTrees$trees[1:10]) #trees_list to run on subset of trees to speed things up for testing
#dev.off()

#perform binary species subset match (SSM) permulation - WITHOUT PLOTTING
permSSM<-getPermsBinaryQuiet(numperms=nperms, fg_vec=fg, sisters_list=sisters_fg, root_sp=myRoot, RERmat=myRER, trees=RERtrees, mastertree=mt, permmode="ssm")

#Produces permulated trees using SSM permulations for one gene tree
#ENSMUSP00000011262 is tree number 611; this "loop" will return an actual "phylo" object so that ape functions will work
#for(i in myTrees$trees[611]){
#	thisTree<-i
#}
#print(thisTree)
#tree_rep = lapply(1:nperms,rep_tree,tree=thisTree)
#root_sp=tree_rep[[1]]$tip.label[1]
##root_sp=sample(tree_rep[[1]]$tip.label, 1)
##gets sisters list for each gene tree, to account for missing data and possible different sisters structure in each gene tree
#temp_sis=gtSistersList(thisTree, fg)
#pdf(paste0("perm_trees.SSMperms",nperms,".one_gene.pdf"), onefile=TRUE, height=15, width=15)
##permulated.binphens = lapply(tree_rep,simBinPhenoSSM,trees=RERtrees,root_sp=root_sp,fg_vec=fg,sisters_list=temp_sis,pathvec=pathvec, plotTreeBool=T)
#permulated.binphens = lapply(tree_rep,simBinPhenoSSM_fromMasterTree,trees=RERtrees,root_sp=root_sp,fg_vec=fg,sisters_list=temp_sis,pathvec=pathvec, plotTreeBool=T)
#dev.off()

#treeSSM<-simBinPhenoSSM(myTrees, masterTree, myRoot, fg, sisters_fg, pathvec, plotTreeBool=F)
#pdf(paste("RERconverge_output.logRTM_binary",dataset,"permSSM_tree.pdf", sep="."))
#plotTreeHighlightBranches(treeSSM,hlspecies=which(treeSSM$edge.length==1), hlcols="red")
#dev.off()


print("Saving permulation output")
save(permSSM, file=paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.parallel",itn,".rds"))
#save(tree_rep, root_sp, temp_sis, permulated.binphens, file=paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM_TEST.one_gene.rds"))
#FOR SISTERS NULL
#save(tree_rep, root_sp, permulated.binphens, file=paste0(myPath, "RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM_TEST.one_gene.rds"))

print("Done with 09_permulations_only.binary.SSM.r")
