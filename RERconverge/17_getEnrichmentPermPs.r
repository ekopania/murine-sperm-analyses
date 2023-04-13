#PURPOSE: Calculate enrichment perm Ps from enrichment perm results

library(topGO)
library(RERconverge)
source("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/GO_permulation.r")

print("Loading data...")
#Load enrichment perm output
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.enrichment.customSisters.rds")
#Load regular permulation output
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.customSisters.rds")
#Load RERconverge run (before permulations)
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")

#Combine permulation and enrichment permulation info
permWenrich<-combpermSSM
permWenrich$enrichStat<-enrich_perms$enrichStat
permWenrich$enrichP<-enrich_perms$enrichP

#Calculate real stats and real enrichments from RERconverge run
print("Calculating real GO enrichment results...")
realstat<-sign(res$Rho)*-log10(res$P)
#Convert rownames to gene IDs for GO enrichment test
prots<-unlist(sapply(rownames(res), function(x) gsub("-.*","",x)))
genes<-protToGene(prots)
names(realstat)<-genes
realenrich<-accGO(na.omit(realstat), 250)
print(dim(realenrich))
print(head(realenrich))

#Get GO permulation P values
print("Calculating permulation P values for GO terms...")
permpGO<-permpvalGOenrich(realenrich, permWenrich, progress=TRUE)
#Correct for multiple tests
permpGO_adj<-p.adjust(permpGO, method="BH")

#Get GO permulation p values in the correct order and append to output file
print("Generating final table...")
permpGO_matched<-permpGO[match(realenrich$GO.ID, names(permpGO))]
permpGO_adj_matched<-permpGO_adj[match(realenrich$GO.ID, names(permpGO_adj))]
stopifnot(all.equal(names(permpGO_matched), realenrich$GO.ID))
stopifnot(all.equal(names(permpGO_adj_matched), realenrich$GO.ID))
stopifnot(all.equal(names(permpGO_adj_matched), names(permpGO_matched)))

GO_with_perm<-as.data.frame(cbind(realenrich, permp.GO=unlist(permpGO_matched), permp.GO.adj=unlist(permpGO_adj_matched)))
print(dim(GO_with_perm))
print(head(GO_with_perm))
write.csv(GO_with_perm, file="RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.customSisters.GOenrichmentPerms.csv", row.names=FALSE)

