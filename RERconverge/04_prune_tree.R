#PURPOSE: Use APE to remove tips that do not have phenotype data

library(ape)

print("Reading in data...")
myTree<-read.tree("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/full_coding_iqtree_astral_rooted_bl.cf.tree")
mySamples<-scan("RTM_taxa_list.txt", what=character(), comment.char="#")
output_treefile<-"pruned_tree.RTM_taxa.tree"
print(myTree)
print(head(mySamples))

print("Pruning...")
outTree<-keep.tip(myTree, mySamples)

print("Pruned tree info:")
print(outTree)

write.tree(outTree, file = output_treefile)

print("Done with 04_prune_tree.R")
