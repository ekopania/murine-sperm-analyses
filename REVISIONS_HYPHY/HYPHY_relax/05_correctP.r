#PURPOSE: Correct P values for multiple tests, sort, get significant loci
mydata<-read.csv("relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.csv", header=TRUE, comment.char="#")
p.adj<-p.adjust(mydata$p.value, method="BH")
newdata<-as.data.frame(cbind(mydata$id, mydata$LRT, mydata$p.value, p.adj, mydata$K))
colnames(newdata)<-c("protID", "LRT", "p.value", "p.adj.BH", "K")
#print(head(newdata))
newdata$LRT<-as.numeric(newdata$LRT)
newdata$p.value<-as.numeric(newdata$p.value)
newdata$p.adj.BH<-as.numeric(newdata$p.adj.BH)
newdata$K<-as.numeric(newdata$K)
sortdat<-newdata[order(newdata$p.adj.BH),]
print(head(sortdat))
print(paste("Total number of loci:", nrow(sortdat)))
write.csv(sortdat, "relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.BHcorrection.csv", row.names=FALSE, quote=FALSE)
sigonly<-sortdat[which(sortdat$p.adj.BH < 0.05),]
print(paste("Number of significant loci:", nrow(sigonly)))
write.csv(sigonly, "relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.BHcorrection.sigOnly.csv", row.names=FALSE, quote=FALSE)

#Get proteins w/ evidence for relaxed purifying selection only
#relaxed<-sortdat[which(sortdat$K < 1),]
#write.csv(relaxed, "relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.BHcorrection.relaxedOnly.csv", row.names=FALSE, quote=FALSE)
relaxed_sig<-sigonly[which(sigonly$K < 1),]
print(paste("Number of significantly relaxed loci:", nrow(relaxed_sig)))
write.csv(relaxed_sig, "relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.BHcorrection.sigOnly.relaxedOnly.csv", row.names=FALSE, quote=FALSE)

#Get proteins w/ intensified selection only
intense_sig<-sigonly[which(sigonly$K > 1),]
print(paste("Number of significantly intensified loci:", nrow(intense_sig)))
write.csv(intense_sig, "relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.BHcorrection.sigOnly.intensifiedOnly.csv", row.names=FALSE, quote=FALSE)
