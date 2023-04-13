#PURPOSE: Use APE to remove tips that do not have phenotype data and make tree ultrametric

library(ape)

print("Reading in data...")
myTree<-read.tree("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/full_coding_iqtree_astral_rooted_bl.cf.tree")
myPhenos<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/murine_repro_pheno_data_ALL.formatted.csv", header=TRUE)
output_treefile<-"pruned_ultrametric_tree.RTM.full_coding_iqtree_astral_rooted_bl.cf.tree"
print(myTree)
print(head(myPhenos))
print(dim(myPhenos))
myPhenos_RTM<-myPhenos[which(myPhenos$RTM!=""),]
print("Species with RTM data only:")
print(dim(myPhenos_RTM))
#print(myPhenos$RTM)

print("Pruning...")
keep_list<-myPhenos_RTM$tree_sample_name
pTree<-keep.tip(myTree, keep_list)

print("Pruned tree info:")
print(pTree)

#Hacky way to deal with no branch lengths leading to tips in ASTRAL output; just changing these to very small branch length
pTree$edge.length<-unlist(sapply(pTree$edge.length, function(x) if(is.na(x)) 0 else x))
pTree$edge.length <- pTree$edge.length + 1e-08 
#print(pTree$edge.length)

outTree<-chronos(pTree)
print("Ultrametric tree info:")
print(outTree)
print("Is tree ultrametric? Should be TRUE")
print(is.ultrametric(outTree))

write.tree(outTree, file = output_treefile)

print("Done with 01_prune_tree.R")
