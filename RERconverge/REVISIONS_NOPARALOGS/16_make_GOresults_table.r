#PURPOSE: Organize a table for GO enrichment results for initial RERconverge run and permulation run

#Read in data
rerc<-read.csv("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.GOenrichment.acceleratedONLY.top500.csv", header=TRUE)
perm<-read.csv("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.permulationsSSM.customSSM.fromMasterTree.GOenrichment.acceleratedONLY.top500.csv", header=TRUE)

#Sort RERconverge df to match permulation df
rerc_sort<-rerc[match(perm$GO.ID, rerc$GO.ID),]
print(head(rerc_sort))
#Check order
stopifnot(all.equal(rerc_sort$GO.ID, perm$GO.ID))

#Get fold-enrichment (obs / exp)
rerc_sort_fe<-rerc_sort$Significant / rerc_sort$Expected
perm_fe<-perm$Significant / perm$Expected

#Set up final dataframe and write to csv
final_df<-as.data.frame(cbind(rerc_sort$GO.ID, rerc_sort$Term, rerc_sort$Annotated, rerc_sort$Significant, rerc_sort$Expected, rerc_sort_fe, rerc_sort$classFisher, rerc_sort$p.adj, perm$Significant, perm$Expected, perm_fe, perm$classFisher, perm$p.adj))
colnames(final_df)<-c("GO.ID", "Term", "Annotated", "Significant", "Expected", "Fold-enrichment", "Pval", "P.adj", "Perm_Significant", "Perm_Expected", "Perm_Fold-enrichment", "Perm_Pval", "Perm_P.adj")

write.csv(final_df, file="GOsummary.RERandPerm.csv", row.names=FALSE, quote=TRUE)
