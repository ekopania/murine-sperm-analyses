#PURPOSE: Correct for multiple tests on positive selection chi2 test

print("Running R script to correct for multiple tests and identify significant genes")

myData<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.bsTest.top10.txt", header=TRUE, col.names=c("protID","lnL","lnL.null","lnL.m1a","LRT","P","LRT_relaxed","P_relaxed"))

myData$p.adj<-p.adjust(myData$P, method="BH")
myData$p_relaxed.adj<-p.adjust(myData$P_relaxed, method="BH")

myData<-myData[, c("protID","lnL","lnL.null","lnL.m1a","LRT","P","p.adj","LRT_relaxed","P_relaxed","p_relaxed.adj")]

write.table(myData, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.bsTest.top10.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

print(paste("Number of significant genes, branch-site test:", length(which(myData$p.adj < 0.05))))
print(paste("Number of significant genes, relaxed selection:", length(which(myData$p_relaxed.adj < 0.05))))

myData_sig<-myData[which(myData$p.adj < 0.05),]

write.table(myData_sig, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.sigOnly.bsTest.top10.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

myData_pos<-myData_sig[which(myData_sig$LRT >= 0),]

write.table(myData_pos, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.sigOnly.posOnly.bsTest.top10.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
