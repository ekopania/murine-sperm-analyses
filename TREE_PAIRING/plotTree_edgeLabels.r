#PURPOSE: Plot pruned tree (species w/ RTM data only) with edge labels to determine which edges to exempt from paring

library(ape)
#myTree<-read.tree("/mnt/beegfs/ek112884/murinae/OU_MODELS/pruned_ultrametric_tree.RTM.full_coding_iqtree_astral.cf.rooted.tree")
myTree<-read.tree("pared_tree.pruned_ultrametric_RTM.keepOUshifts.bl100.tree/iter-1-pared.tre")
print(myTree)
#pdf("pruned_ultrametric_RTM.withLabels.pdf", height=25, width=25)
pdf("pared_tree.pruned_ultrametric_RTM.keepOUshifts.bl100.tree/tree_plot.iter-1-pared.pdf", height=25, width=25)
plot(myTree)
nodelabels()
dev.off()
