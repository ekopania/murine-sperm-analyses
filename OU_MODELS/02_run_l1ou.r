#PURPOSE: Run l1ou to identify convergent shifts in optima for a trait across the murine tree

library(l1ou)

print("Reading in data...")
myTree<-read.tree("pruned_ultrametric_tree.RTM.full_coding_iqtree_astral_rooted_bl.cf.tree")
myTree$tip.label<-unname(unlist(sapply(myTree$tip.label, function(x) paste(strsplit(x, "_")[[1]][1], strsplit(x,"_")[[1]][2],sep="_"))))
#names(myTree$tip.label)<-c(1:length(myTree$tip.label))
myPhenos<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/murine_repro_pheno_data_RTM.formatted.csv")
myPhenos$short_name<-unlist(sapply(myPhenos$tree_sample_name, function(x) paste(strsplit(x, "_")[[1]][1], strsplit(x, "_")[[1]][2], sep="_")))
normTree<-normalize_tree(myTree, check.ultrametric=TRUE)
print(normTree)
print(head(myPhenos))

print("Adjust data...")
#this_pheno<-myPhenos$RTM #Percent of body weight
this_pheno<-log(myPhenos$RTM) #log(percent of body weight)
pdf("RTM.log_transformed.hist.pdf")
hist(this_pheno, main="log(RTM)")
dev.off()
names(this_pheno)<-myPhenos$short_name
#this_pheno_ordered<-this_pheno[match(names(this_pheno), normTree$tip.label)]
#normTree$tip.label<-unlist(sapply(normTree$tip.label, function(x) paste(strsplit(x, "_")[[1]][1], strsplit(x, "_")[[1]][2], sep="_")))
print(this_pheno)
print(normTree)
#print(normTree$tip.label)
#print(names(this_pheno) %in% normTree$tip.label)
myAdj_data<-adjust_data(normTree, this_pheno)

print("Estimate shift configuration...")
eModel<-estimate_shift_configuration(myAdj_data$tree, myAdj_data$Y, criterion="AICc", root.model="OUrandomRoot")
print(summary(eModel))

print("Running OU fit...")
myOU<-fit_OU(eModel$tree, eModel$Y, eModel$profile$configurations[[1]], criterion="AICc", root.model="OUrandomRoot")
pdf("OUmodel.RTM.pdf", height=24, width=7.5)
plot(myOU)
dev.off()

print("Testing for convergent regimes...")
fit_conv <- estimate_convergent_regimes(eModel, criterion="AIC", method="rr")
print(fit_conv)
pdf("OUmodel_convergence.RTM.pdf", height=24, width=7.5)
plot(fit_conv)
dev.off()

print("Done with 02_run_l1ou.r")
