#PURPOSE: Combine permulation runs and analyze them

library(RERconverge)

myPath<-"/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/PERMULATIONS_OUTPUT_OUmodel_customSSM_fromMasterTree/"

print("Loading data...")
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")

perm.files<-list.files(path=myPath, pattern=".permulationsSSM.parallel\\d+.rds")
print("Permulation files to combine:")
print(perm.files)

combpermSSM<-c()
for(f in perm.files){
	num<-sub("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.parallel", "", sub(".rds", "", f))
	print(num)
	load(paste0(myPath,f))
	assign(paste0("permSSM", num), permSSM)
	if(length(combpermSSM)==0){
		combpermSSM<-permSSM
	} else{
		combpermSSM<-combinePermData(combpermSSM, permSSM, enrich=F)
	}
}
#warnings()

#print("Testing for  perms:")
#list.<-sapply(c(1:), function(x) paste0("permSSM",x))
#print(list.)
#print("Combining first  permulations:")
#combpermSSM<-combinePermData(list(list.), enrich=F)
print("Getting pvals")
permpvalSSM<-permpvalcor(res,combpermSSM)
print(head(permpvalSSM))
print(length(which(permpvalSSM$permpval < 0.05)))
print(median(permpvalSSM$permpval, na.rm=TRUE))

write.table(permpvalSSM, file="RERconverge_output.logRTM_binary_permpvalSSM.customSSM.OUmodel.RTMspeciesOnly.fromMasterTree.txt")

# add permulations to real results
print("Appending p-vals")
print(head(res))
print(head(permpvalSSM))
res$permpval<-permpvalSSM$permpval[match(rownames(res), rownames(permpvalSSM))]
res$permpvaladj<-p.adjust(res$permpval, method="BH")
res$permstats<-permpvalSSM$permstats[match(rownames(res), rownames(permpvalSSM))]
write.csv(res, file="RERconverge_output.logRTM_binary_residuals.OUmodel.RTMspeciesOnly.withPermulationsSSM.customSSM.fromMasterTree.csv", row.names=TRUE)

#Save data
save(combpermSSM, permpvalSSM, file="RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.fromMasterTree.rds")
