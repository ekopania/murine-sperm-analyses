#PURPOSE: Run RERconverge correlations on log-transformed phenotype data

library(RERconverge)
print("Correlating RERs with log-transformed phenotype data")

print("Loading data...")
#Load data output from RERconverge
#We will keep the RER values since these are independent of phenotype, but replace charpaths and res using log-transformed pheno data
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")
print("Loaded variables:")
ls()
#print(myTrees)
#Read in phenotype data
myData<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/murine_repro_pheno_data_RTM.formatted.csv", header=TRUE)
myPheno<-log(myData$RTM)
names(myPheno)<-myData$tree_sample_name

print("Generating trait tree...")
#generate trait tree
charpaths<-char2Paths(myPheno, myTrees)

print("Correlate...")
#calculate correlation statistics
res<-correlateWithContinuousPhenotype(myRER, charpaths, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
print(head(res))
#order correlation statistics by corrected p-value and save in csv file
res_pSort<-res[order(res$p.adj),]
write.csv(res_pSort, file="RERconverge_output.logRTM_continuous_residuals.RTMspeciesOnly.csv", row.names=TRUE)

#Save data
save(myRER, charpaths, res, file="RERconverge_output.logRTM_continuous.RTMspeciesOnly.rds")

print("Done with 07.1_get_RERs_logPheno.r")
