#PURPOSE: Plot how often certain nodes and tips end up in the foreground in permulated trees

library(RERconverge)

#First 100 trees
#subTrees<-scan("submit/temp.from_masterTree.midpoint_root.10loci.tre", what=character(), nmax=100)
#Next 100 trees; increase skip by 100 for each set of 100 trees
subTrees<-scan("submit/temp.from_masterTree.midpoint_root.10loci.tre", what=character(), skip=900, nmax=100)

fg_nodes<-c()
fg_tips<-c()
for(t in subTrees){
	tree<-read.tree(text=t)
	#print(tree)
	edge = tree$edge
	#print(edge)
	edge.length=tree$edge.length
	#print(edge.length)
	ind.fg.edge = which(edge.length == 1)
	#print(ind.fg.edge)
	nodeIds.fg.edge = edge[ind.fg.edge,] # all foreground edges in the temp tree
	#print(edge[ind.fg.edge,])
	print(nodeIds.fg.edge)
	fg.tips = tree$tip.label[nodeIds.fg.edge[,2]]
	print(fg.tips)
	fg_nodes<-c(fg_nodes, nodeIds.fg.edge[,2])
	fg_tips<-c(fg_tips, fg.tips)
}

uniq_nodes<-unique(fg_nodes)
node_counts<-c()
for(i in uniq_nodes){
        count<-length(which(fg_nodes==i))
        node_counts<-c(node_counts, count)
}
names(node_counts)<-uniq_nodes


uniq_tips<-na.omit(unique(fg_tips))
print(uniq_tips)
tip_counts<-c()
for(i in uniq_tips){
	count<-length(which(fg_tips==i))
	tip_counts<-c(tip_counts, count)
}
names(tip_counts)<-uniq_tips
#Add zeros
zero_tips<-tree$tip.label[which(!(tree$tip.label %in% names(tip_counts)))]
zvec<-rep(0, length(zero_tips))
names(zvec)<-zero_tips
tip_counts<-c(tip_counts, zvec)
print(tip_counts)

print(paste("Number of unique nodes in the foreground in at least one permulated tree:", length(node_counts)))
print(paste("Number of unique tips in the foreground in at least one permulated tree:", sum(tip_counts!=0)))

pdf("perm_foregrounds.pdf", onefile=TRUE, height=12, width=12)
barplot(node_counts[sort(names(node_counts))], main="Foreground nodes", xlab="Node labels", ylab="Number of permulated trees", las=2)
par(mar = c(20.1, 4.1, 4.1, 2.1))
barplot(tip_counts[sort(names(tip_counts))], main="Foreground tips", ylab="Number of permulated trees", las=2)
dev.off()
