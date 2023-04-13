#PURPOSE: Correct for multiple tests on positive selection chi2 test

print("Running R script to correct for multiple tests and identify significant genes")

myData<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.m1aVm2a.txt", header=FALSE, col.names=c("protID","lnL","lnL.null","diff","P"))

myData$p.adj<-p.adjust(myData$P, method="BH")

write.table(myData, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.pValCor.m1aVm2a.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

print(paste("Number of significant genes:", length(which(myData$p.adj < 0.05))))

myData_sig<-myData[which(myData$p.adj < 0.05),]

write.table(myData_sig, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.pValCor.sigOnly.m1aVm2a.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

myData_pos<-myData_sig[which(myData_sig$diff >= 0),]

write.table(myData_pos, file="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.pValCor.sigOnly.posOnly.m1aVm2a.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
