#PURPOSE: Run gls function in nlme package to test for significant associations among reproductive phenotypes while controling for phylogeny
#Based on this tutorial http://www.phytools.org/Cordoba2017/ex/4/PGLS.html to test out nmle package
#NOTE: Because of an update in nlme, the form parameter must be specified when setting correlation (in contrast to the tutorial); form is the "taxa covariate" and should usually be a list of species names that correspond to those in the tree

library(ape)
library(nlme)
library(phytools)
library(phylolm)
library(ggplot2)

#Read in data
#pahl<-read.csv("murine_repro_pheno_data_PahlEtAl2018.formatted.csv",header=TRUE)
#mcl<-read.csv("murine_repro_pheno_data_McLennanEtAl2016.formatted.csv",header=TRUE)
rtm<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/murine_repro_pheno_data_RTM.formatted.csv",header=TRUE)
spermMorhoAll<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/murine_repro_pheno_data_ALL.formatted.csv",header=TRUE)
#print(colnames(pahl))
#print(colnames(mcl))
mycols<-c("tree_sample_name", "Head_length_um", "Head_width_um", "Head_area_um.2", "Apical_hook_length_um", "Apical_hook_angle", "Ventral_process_length_um", "Ventral_process_angle", "Midpiece_length_um", "Principal_and_end_piece_length_um", "Tail_length_um", "num_hooks_or_vp")
#NOTE: Need to fill in num_hooks_or_vp later
#pahl_phenos<-pahl[, which(colnames(pahl) %in% mycols)]
#mcl_phenos<-mcl[, which(colnames(mcl) %in% mycols)]
spm_phenos<-spermMorhoAll[, which(colnames(spermMorhoAll) %in% mycols)]
#print(colnames(pahl_phenos))
#print(colnames(mcl_phenos))

#pahl_all<-as.data.frame(cbind(pahl_phenos, Head_width_um=rep(NA,nrow(pahl_phenos)), Ventral_process_length_um=rep(NA,nrow(pahl_phenos)), Ventral_process_angle=rep(NA,nrow(pahl_phenos)), num_hooks_or_vp=rep(NA,nrow(pahl_phenos))))
#mcl_all<-as.data.frame(cbind(mcl_phenos, Head_area_um.2=rep(NA,nrow(mcl_phenos))))
spm_all<-as.data.frame(spm_phenos)
#print(dim(pahl_all))
#print(dim(mcl_all))
#print(colnames(pahl_all))
#print(colnames(mcl_all))
#pahl_final<-pahl_all[, order(match(colnames(pahl_all), mycols))]
#mcl_final<-mcl_all[, order(match(colnames(mcl_all), mycols))]
spm_final<-spm_all[, order(match(colnames(spm_all), mycols))]
#print(head(pahl_final))
#print(head(mcl_final))
#sperm_pheno<-as.data.frame(rbind(pahl_final, mcl_final))

#Add hook presence/absence as binary trait
print(head(spm_final))
hook_pa<-sapply(spm_final$num_hooks_or_vp, function(x) if(is.na(x)) NA else if(x>0) 1 else 0)

#final DF
sperm_pheno<-as.data.frame(cbind(spm_final, hook_pa))
print(dim(sperm_pheno))
print(head(sperm_pheno))
print(tail(sperm_pheno))

mytree<-read.tree("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/full_coding_iqtree_astral_rooted_bl.cf.tree")
#mytree<-read.tree("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/OU_MODELS/pruned_ultrametric_tree.RTM.full_coding_iqtree_astral_rooted_bl.cf.tree")
#ASTRAL output does not give branch lengths on terminal branches, so R reads them in as "NaN"; replace with "1"
new_edges<-as.numeric(sapply(mytree$edge.length, function(x) sub(NaN,1,x)))
mytree$edge.length<-new_edges

#Make sure data and tree have same taxa with same names
keep<-Reduce(intersect, list(sperm_pheno$tree_sample_name, rtm$tree_sample_name, mytree$tip.label))
print(length(keep))
keepSperm<-sperm_pheno[which(sperm_pheno$tree_sample_name %in% keep), ]
keepRTM<-rtm[which(rtm$tree_sample_name %in% keep), ]
keepSperm$tree_sample_name<-as.character(keepSperm$tree_sample_name)
keepRTM$tree_sample_name<-as.character(keepRTM$tree_sample_name)
keepRTM_sorted<-keepRTM[order(match(keepRTM$tree_sample_name, keepSperm$tree_sample_name)), ]
print(dim(keepSperm))
#print(keepSperm)
print(dim(keepRTM_sorted))
stopifnot(all.equal(keepSperm$tree_sample_name, keepRTM_sorted$tree_sample_name))
mydata<-as.data.frame(cbind(keepSperm, rtm=keepRTM_sorted$RTM))
print(mydata)

#Loop through every combination of traits
#"Head_length_um" and "Midpiece_length_um" not converging with some traits, so removing for now
#mytraits<-c("rtm", "Head_length_um", "Head_width_um", "Head_area_um.2", "Apical_hook_length_um", "Apical_hook_angle", "Ventral_process_length_um", "Ventral_process_angle", "Principal_and_end_piece_length_um", "Midpiece_length_um", "Tail_length_um", "hook_pa")
#mytraits<-c("rtm", "Head_width_um", "Head_area_um.2", "Apical_hook_length_um", "Apical_hook_angle", "Ventral_process_length_um", "Ventral_process_angle", "Principal_and_end_piece_length_um", "Tail_length_um", "hook_pa")
mytraits<-c("rtm", "Apical_hook_length_um", "Apical_hook_angle", "hook_pa")
physig<-c()
physig_ps<-c()
models<-c()
bm_ps<-c()
pagel_ps<-c()
lambdas<-c()
ou_ps<-c()
samp_sizes<-c()
pdf("trait_regressions.logPlus1Transform.pdf", onefile=TRUE)
#pdf("trait_regressions.pdf", onefile=TRUE)
for(i in 1:(length(mytraits)-1)){
	#Get phylogenetic signal for trait
	temp_pheno<-mydata[, mytraits[i]]
	names(temp_pheno)<-mydata$tree_sample_name
	print(paste("Phylogenetic signal for trait", mytraits[i]))
	print(paste("Number of taxa with trait data:", length(which(!(is.na(temp_pheno))))))
	print(head(temp_pheno))
	ps<-phylosig(mytree, temp_pheno, method="lambda", test=TRUE)
	print(ps)
	physig_ps<-c(physig_ps, ps$P)
	physig<-c(physig, ps$lambda)
	for(j in (i+1):length(mytraits)){
		t1<-mytraits[i]
		t2<-mytraits[j]
		print(paste("Running model", t1, "~", t2))
		#Remove taxa with NA for the traits of interest; commented sections also remove zeros
		if(t1=="hook_pa"){
			#t1_keep<-intersect(which(!(is.na(mydata[,t1]))), which(mydata[,t1] != 0))
			t2_keep<-intersect(which(!(is.na(mydata[,t2]))), which(mydata[,t2] != 0))
			t1_keep<-which(!(is.na(mydata[,t1])))
			#t2_keep<-which(!(is.na(mydata[,t2])))
		} else if(t2=="hook_pa"){
			t1_keep<-intersect(which(!(is.na(mydata[,t1]))), which(mydata[,t1] != 0))
			#t2_keep<-intersect(which(!(is.na(mydata[,t2]))), which(mydata[,t2] != 0))
			#t1_keep<-which(!(is.na(mydata[,t1])))
			t2_keep<-which(!(is.na(mydata[,t2])))
		} else{
			t1_keep<-intersect(which(!(is.na(mydata[,t1]))), which(mydata[,t1] != 0))
                        t2_keep<-intersect(which(!(is.na(mydata[,t2]))), which(mydata[,t2] != 0))
			#t1_keep<-which(!(is.na(mydata[,t1])))
			#t2_keep<-which(!(is.na(mydata[,t2])))
		}
		tmp_data<-mydata[intersect(t1_keep, t2_keep), c("tree_sample_name",t1,t2)]
		if(nrow(tmp_data)<=1){
			print("One or no samples with data for both these traits; skipping...")
			next #Go to next pair of traits if there are no samples with overlapping data
		}
		pruned.tree<-drop.tip(mytree, mytree$tip.label[which(!(mytree$tip.label %in% tmp_data$tree_sample_name))])
		stopifnot(all.equal(sort(pruned.tree$tip.label), sort(tmp_data$tree_sample_name)))

		print(tmp_data)
		print(pruned.tree)
		
		#Log-transform pheno data following Pahl et al. 2018; add 1 for traits with value 0
		log_transform<-function(x){
			log(as.numeric(as.character(x)) + 1)
		}
		log_data<-mapply(log_transform, tmp_data[,2:3])
		merge_data<-as.data.frame(cbind(tmp_data$tree_sample_name, log_data))
		#merge_data<-as.data.frame(cbind(tmp_data$tree_sample_name, tmp_data[,2:3]))
		colnames(merge_data)<-c("tree_sample_name", "trait1", "trait2")
		final_data<-merge_data[order(match(merge_data$tree_sample_name, pruned.tree$tip.label)),]
		final_data$trait1<-as.numeric(as.character(final_data$trait1))
		final_data$trait2<-as.numeric(as.character(final_data$trait2))
		n<-nrow(final_data)
		print(final_data)
		#Don't try to run lm if there is no variance
		if( (length(unique(final_data$trait1))==1) || (length(unique(final_data$trait2))==1) ){
			next
		} else{
			#add to sample size vector
			samp_sizes<-c(samp_sizes, n)
			#Run models
	
			#Simple BM model to test for relationship between two phenotypic variables (Lnote and Lnalt, song note length and altitude)
			#First variable in corBrownian() is the gamma parameter; 1 is the default
			#Newest version of nlme requires the "form" parameter to run properly; should be a list of taxa in the tree
			fm<-as.formula("trait1~trait2")
			print(fm)
			models<-c(models, paste0(t1,"~",t2))
			bm_model<-gls(fm, data=final_data, correlation=corBrownian(1, pruned.tree, form=~tree_sample_name))
			print(summary(bm_model))
			bm_tab<-as.data.frame(summary(bm_model)$tTable)
			bm_ps<-c(bm_ps, bm_tab$p[2])
			
			#Repeat w/ Pagel's lambda model instead of BM model 
			#Notice the "form" parameter is added to all the corPagel() functions below
			pa_model<-gls(fm, data=final_data, correlation=corPagel(1, pruned.tree, form=~tree_sample_name))
			print(summary(pa_model))
			pa_tab<-as.data.frame(summary(pa_model)$tTable)
			pagel_ps<-c(pagel_ps, pa_tab$p[2])
			pa_cor<-summary(pa_model)$corBeta
			print(pa_cor)
			lambdas<-c(lambdas, pa_model$modelStruct$corStruct) 
	
			#Repeat w/ OU model
			#Scaling branch lengths by a constant (*100) otherwise model doesn't converge
			#See https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/
			tempTree<-pruned.tree
			tempTree$edge.length<-tempTree$edge.length * 100
			ou_model<-gls(fm, data=final_data, correlation=corMartins(1, tempTree, form=~tree_sample_name))
			print(summary(ou_model))
			ou_tab<-as.data.frame(summary(ou_model)$tTable)
                        ou_ps<-c(ou_ps, ou_tab$p[2])
			
			#Log likelihood BM vs Pagel
			
			#Logistic regression (Garland and Ives 2010) - better for discrete traits (# hooks)
			if( (t1=="hook_pa") || (t2=="hook_pa") ){
				print("Logistic regression for binary trait:")
				#Response variable must be discrete (num_hooks_or_vp has to be trait1)
				if(t1=="hook_pa"){
					phyloglm_df<-as.data.frame(tmp_data[,2:3])
				} else{
					phyloglm_df<-as.data.frame(cbind(tmp_data[,3], tmp_data[,2]))
				}
				rownames(phyloglm_df)<-tmp_data$tree_sample_name
				colnames(phyloglm_df)<-c("trait1","trait2")
				print(head(phyloglm_df))
				#btol = 20 to address warning
				gi_model<-phyloglm(fm, data=phyloglm_df, phy=pruned.tree, method="logistic_IG10", btol=30)
				print(summary(gi_model))
				#Plot
				p<-ggplot(phyloglm_df, aes(x=trait2, y=trait1)) + geom_point() + xlim(-1,5) + theme_minimal()
                        	#Logistic regression
			        p<-p + stat_smooth(method = "glm", col = "red", method.args=list(family=binomial), fullrange=TRUE)
                                if(t1=="hook_pa"){
					p<-p + labs(title=paste(t1, "vs", t2), x=t2, y=t1)
				} else{
					p<-p + labs(title=paste(t2, "vs", t1), x=t1, y=t2)
				}
				#p<-p + theme_minimal() + coord_flip()
				print(p)
			} else{ #Linear regression
				#Plot
				p<-ggplot(final_data, aes(x=trait1, y=trait2)) + geom_point()
				p<-p + stat_smooth(method = "glm", col = "red")
				p<-p + labs(title=paste(t1, "vs", t2), x=paste0("log(",t1,")"), y=paste0("log(",t2,")"))
				p<-p + theme_minimal()
                        	print(p)
			}
		}
	}
}
#dev.off()

#Get phylogenetic signal for last trait (first for loop doesn't get to)
temp_pheno<-mydata[, mytraits[length(mytraits)]]
names(temp_pheno)<-mydata$tree_sample_name
print(paste("Phylogenetic signal for trait", mytraits[length(mytraits)]))
print(paste("Number of taxa with trait data:", length(which(!(is.na(temp_pheno))))))
print(head(temp_pheno))
ps<-phylosig(mytree, temp_pheno, method="lambda", test=TRUE)
physig_ps<-c(physig_ps, ps$P)
physig<-c(physig, ps$lambda)

#Get results for RTM vs head length, because head length was removed from the loop due to convergence issues with some other phenotypes
#t1<-"rtm"
#t2<-"Head_length_um"
#print(paste("Running model", t1, "~", t2))
#t1_keep<-which(!(is.na(mydata[,t1])))
#t2_keep<-which(!(is.na(mydata[,t2])))
#tmp_data<-mydata[intersect(t1_keep, t2_keep), c("tree_sample_name",t1,t2)]
#pruned.tree<-drop.tip(mytree, mytree$tip.label[which(!(mytree$tip.label %in% tmp_data$tree_sample_name))])
#stopifnot(all.equal(sort(pruned.tree$tip.label), sort(tmp_data$tree_sample_name)))
#print(tmp_data)
#print(pruned.tree)
##Log-transform pheno data following Pahl et al. 2018; add 1 for traits with value 0
#log_transform<-function(x){
#	log(as.numeric(as.character(x)))# + 1)
#}
#log_data<-mapply(log_transform, tmp_data[,2:3])
#merge_data<-as.data.frame(cbind(tmp_data$tree_sample_name, log_data))
#colnames(merge_data)<-c("tree_sample_name", "trait1", "trait2")
#final_data<-merge_data[order(match(merge_data$tree_sample_name, pruned.tree$tip.label)),]
#final_data$trait1<-as.numeric(as.character(final_data$trait1))
#final_data$trait2<-as.numeric(as.character(final_data$trait2))
#n<-nrow(final_data)
#print(final_data)
##Don't try to run lm if there is no variance
#if( (length(unique(final_data$trait1))==1) || (length(unique(final_data$trait2))==1) ){
#	next
#} else{
#	#add to sample size vector
#        samp_sizes<-c(samp_sizes, n)
#        #Run models
#	#Simple BM model to test for relationship between two phenotypic variables (Lnote and Lnalt, song note length and altitude)
#        fm<-as.formula("trait1~trait2")
#        print(fm)
#        models<-c(models, paste0(t1,"~",t2))
#        bm_model<-gls(fm, data=final_data, correlation=corBrownian(1, pruned.tree, form=~tree_sample_name))
#        print(summary(bm_model))
#        bm_tab<-as.data.frame(summary(bm_model)$tTable)
#        bm_ps<-c(bm_ps, bm_tab$p[2])
#	#Repeat w/ Pagel's lambda model instead of BM model
#        pa_model<-gls(fm, data=final_data, correlation=corPagel(1, pruned.tree, form=~tree_sample_name))
#        print(summary(pa_model))
#        pa_tab<-as.data.frame(summary(pa_model)$tTable)
#        pagel_ps<-c(pagel_ps, pa_tab$p[2])
#        pa_cor<-summary(pa_model)$corBeta
#        print(pa_cor)
#        lambdas<-c(lambdas, pa_model$modelStruct$corStruct)
#
#        #Repeat w/ OU model
#        tempTree<-pruned.tree
#        tempTree$edge.length<-tempTree$edge.length * 100
#        ou_model<-gls(fm, data=final_data, correlation=corMartins(1, tempTree, form=~tree_sample_name))
#        print(summary(ou_model))
#        ou_tab<-as.data.frame(summary(ou_model)$tTable)
#        ou_ps<-c(ou_ps, ou_tab$p[2])
#
#	#Plot
#	p<-ggplot(final_data, aes(x=trait1, y=trait2)) + geom_point()
#	p<-p + stat_smooth(method = "glm", col = "red")
#	p<-p + labs(title=paste(t1, "vs", t2), x=paste0("log(",t1,")"), y=paste0("log(",t2,")"))
#	p<-p + theme_minimal()
#	print(p)
#}
dev.off()

#Collect p values, p.adjust, save in nice csv or tsv format
bm_adj<-p.adjust(bm_ps, method="fdr")
pagel_adj<-p.adjust(pagel_ps, method="fdr")
ou_adj<-p.adjust(ou_ps, method="fdr")
#mystats<-as.data.frame(cbind(models, bm_ps, bm_adj, pagel_ps, pagel_adj))
mystats<-as.data.frame(cbind(models, bm_ps, bm_adj, pagel_ps, pagel_adj, lambdas, ou_ps, ou_adj, samp_sizes))
#print(mystats)
write.table(mystats, file="gls_p_values.logPlus1Transform.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#write.table(mystats, file="gls_p_values.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

physig_adj<-p.adjust(physig_ps, method="fdr")
physig_stats<-as.data.frame(cbind(mytraits, lambda=physig, physig_ps, physig_adj))
write.table(physig_stats, file="phylogenetic_signal.logPlus1Transform.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#write.table(physig_stats, file="phylogenetic_signal.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
