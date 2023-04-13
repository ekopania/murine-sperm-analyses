#PURPOSE: Use phytools to plot trait divergence across a tree

library(phytools)

#Read in data
pahl<-read.csv("murine_repro_pheno_data_PahlEtAl2018.formatted.csv",header=TRUE)
mcl<-read.csv("murine_repro_pheno_data_McLennanEtAl2016.formatted.csv",header=TRUE)
rtm<-read.csv("murine_repro_pheno_data_RTM.formatted.csv",header=TRUE)
spm<-read.csv("murine_repro_pheno_data.pheno_detail_set.csv",header=TRUE)
mycols<-c("tree_sample_name", "Head_length_um", "Head_width_um", "Head_area_um.2", "Apical_hook_length_um", "Apical_hook_angle", "Ventral_process_length_um", "Ventral_process_angle", "Midpiece_length_um", "Principal_and_end_piece_length_um", "Tail_length_um", "num_hooks_or_vp")
#NOTE: Need to fill in num_hooks_or_vp later
pahl_phenos<-pahl[, which(colnames(pahl) %in% mycols)]
mcl_phenos<-mcl[, which(colnames(mcl) %in% mycols)]
spm_phenos<-as.data.frame(spm[, which(colnames(spm) %in% mycols)])

pahl_all<-as.data.frame(cbind(pahl_phenos, Head_width_um=rep(NA,nrow(pahl_phenos)), Ventral_process_length_um=rep(NA,nrow(pahl_phenos)), Ventral_process_angle=rep(NA,nrow(pahl_phenos)), num_hooks_or_vp=rep(NA,nrow(pahl_phenos))))
mcl_all<-as.data.frame(cbind(mcl_phenos, Head_area_um.2=rep(NA,nrow(mcl_phenos))))
pahl_final<-pahl_all[, order(match(colnames(pahl_all), mycols))]
mcl_final<-mcl_all[, order(match(colnames(mcl_all), mycols))]
spm_final<-spm_phenos[, order(match(colnames(spm_phenos), mycols))]
#sperm_pheno<-as.data.frame(rbind(pahl_final, mcl_final))
sperm_pheno<-spm_final
print(dim(sperm_pheno))
print(head(sperm_pheno))

mytree<-read.tree("/mnt/beegfs/ek112884/murinae/murine_species_tree.noRef.astralMP_output.rooted.tre")
#ASTRAL output does not give branch lengths on terminal branches, so R reads them in as "NaN"; replace with "1"
new_edges<-as.numeric(sapply(mytree$edge.length, function(x) sub(NaN,1,x)))
mytree$edge.length<-new_edges

#Make sure data and tree have same taxa with same names
#keep<-Reduce(intersect, list(sperm_pheno$tree_sample_name, rtm$tree_sample_name, mytree$tip.label))
#keepSperm<-sperm_pheno[which(sperm_pheno$tree_sample_name %in% keep), ]
#keepRTM<-rtm[which(rtm$tree_sample_name %in% keep), ]
#keepSperm$tree_sample_name<-as.character(keepSperm$tree_sample_name)
#keepRTM$tree_sample_name<-as.character(keepRTM$tree_sample_name)
#keepRTM_sorted<-keepRTM[order(match(keepRTM$tree_sample_name, keepSperm$tree_sample_name)), ]
#print(dim(keepSperm))
#print(dim(keepRTM))
#stopifnot(all.equal(keepSperm$tree_sample_name, keepRTM_sorted$tree_sample_name))
#mydata<-as.data.frame(cbind(keepSperm, rtm=keepRTM_sorted$RTM))

#File containing all trees
pdf("trait_trees.phytools.pdf", onefile=TRUE, height=11, width=8.5)
#Loop through every trait
#These traits cause error for some reason; need to figure out why: Apical_hook_angle", "Ventral_process_angle", "Midpiece_length_um", "Principal_and_end_piece_length_um", "Tail_length_um"
mytraits<-c("Head_length_um", "Head_width_um", "Head_area_um.2", "Apical_hook_length_um", "Ventral_process_length_um", "num_hooks_or_vp", "RTM")
for(t in mytraits){
	print(paste("Working on trait", t))
	#Remove taxa with NA or 0 for the traits of interest
	if(t=="RTM"){
		keep<-intersect(which(!(is.na(rtm[,t]))), which(rtm[,t] != 0))
		tmp_data<-rtm[keep, c("tree_sample_name",t)]
	} else{
		keep<-intersect(which(!(is.na(sperm_pheno[,t]))), which(sperm_pheno[,t] != 0))
		tmp_data<-sperm_pheno[keep, c("tree_sample_name",t)]
	}
	if(nrow(tmp_data)==0){
		next #Go to next trait if there are no samples with data
	}
	pruned.tree<-drop.tip(mytree, mytree$tip.label[which(!(mytree$tip.label %in% tmp_data$tree_sample_name))])
	short_tips<-sapply(pruned.tree$tip.label, function(x) gsub("_.*$","",gsub("-.*$","",x)))
	pruned.tree$tip.label<-short_tips
	filt_data<-tmp_data[which(tmp_data$tree_sample_name %in% mytree$tip.label),]
	short_name<-sapply(filt_data$tree_sample_name, function(x) gsub("_.*$","",gsub("-.*$","",x)))
	stopifnot(all.equal(sort(pruned.tree$tip.label), sort(short_name)))
	
	vec_data<-as.numeric(as.character(filt_data[, t]))
	names(vec_data)<-short_name
	print(vec_data)
	print(pruned.tree)
	
	mymin<-min(vec_data)-1
	mymax<-max(vec_data)+1
	print(vec_data)
	print(mymin)
	print(mymax)
	if(t=="Head_area_um.2"){ #Ancestral state CIs extra large for this one so just hard coding
		obj<-contMap(pruned.tree, vec_data, plot=FALSE, lims=c(6, 45))
	} else{
		obj<-contMap(pruned.tree, vec_data, plot=FALSE, lims=c(mymin,mymax))
	}
	obj<-setMap(obj,colors=c("blue","orange"))
	plot(obj, lwd=7, xlim=c(-0.5,30))
	title(main=paste("Tree for trait:", t))
	errorbar.contMap(obj)
}

dev.off()
