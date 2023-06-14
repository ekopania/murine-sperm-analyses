#PURPOSE: Get stats such as ASTRAL support, gene concordance, and site concordance from the node labels on the tree

library(ape)
#Pruned 78 species repro tree
mytree<-read.tree("pruned_ultrametric_tree.RTM.full_coding_iqtree_astral_rooted_bl.cf.tree")

#Full coding tree
#mytree<-read.tree("../full_coding_iqtree_astral_rooted_bl.cf.tree")

print(mytree)

mylabs<-mytree$node.label
mylabs_split<-sapply(mylabs, function(x) strsplit(x, split="/"))
print("Node labels:")
print(mylabs_split)

print("ASTRAL support values:")
#Get astral support values (first value)
astral_sup<-c()
for(i in 1:length(mylabs_split)){
	astral_sup<-c(astral_sup, mylabs_split[[i]][1])
}
print(astral_sup)

print("Second value - GCF:")
#Get second value - GCF
second_val<-c()
for(i in 1:length(mylabs_split)){
	second_val<-c(second_val, mylabs_split[[i]][2])
}
print(second_val)
print(median(as.numeric(second_val), na.rm=TRUE))

print("Third value - SCF:")
#Get second value - SCF
third_val<-c()
for(i in 1:length(mylabs_split)){
	third_val<-c(third_val, mylabs_split[[i]][3])
}
print(third_val)
print(median(as.numeric(third_val), na.rm=TRUE))

#print("Fourth value:")
#Get fourth value - ???
#fourth_val<-c()
#for(i in 1:length(mylabs_split)){
#        fourth_val<-c(fourth_val, mylabs_split[[i]][4])
#}
#print(fourth_val)
#print(median(as.numeric(fourth_val), na.rm=TRUE))

#print("Fifth value:")
#Get fifth value - ???
#fifth_val<-c()
#for(i in 1:length(mylabs_split)){
 #       fifth_val<-c(fifth_val, mylabs_split[[i]][5])
#}
#print(fifth_val)
#print(median(as.numeric(fifth_val), na.rm=TRUE))
