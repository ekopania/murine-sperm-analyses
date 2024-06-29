#PURPOSE: Get median dN/dS for tissue spec genes, ovary spec genes, compare to testis-spec

ligenes<-read.csv("LiEtal2017_mouse_tissue_specific.withProtID.csv", header=TRUE, comment.char="#")
dnds<-read.table("/projects/ek112884e/murines_revisions/PAML_RERUNS/reproPhenoSet_fromFullCoding_M0_dnds.txt", header=TRUE)

prots0<-ligenes$protID
prots<-unlist(sapply(prots0, function(x) unlist(strsplit(x, split=";"))))

print(head(prots))
print(length(prots))

#All genes versus tissue-specific
ts_dnds<-dnds[which(dnds$protID %in% prots),]
print(paste("Number of tissue-specific proteins in dataset:", nrow(ts_dnds)))
print(paste("Median dN/dS for all tissue-specific proteins:", median(ts_dnds$dN.dS, na.rm=TRUE)))
print("Wilcox test compared to all genes:")
result<-wilcox.test(dnds$dN.dS, ts_dnds$dN.dS)
print(result)

#Ovary specific
ovary_prots0<-ligenes$protID[which(ligenes$tissue=="Ov")]
ovary_prots<-unlist(sapply(ovary_prots0, function(x) unlist(strsplit(x, split=";"))))

print(head(ovary_prots))
print(length(ovary_prots))

ov_dnds<-dnds[which(dnds$protID %in% ovary_prots),]
print(paste("Number of tissue-specific proteins in ovaries:", nrow(ov_dnds)))
print(paste("Median dN/dS for ovary-specific proteins:", median(ov_dnds$dN.dS, na.rm=TRUE)))
print("Wilcox test compared to all genes:")
result<-wilcox.test(dnds$dN.dS, ov_dnds$dN.dS)
print(result)
print("Wilcox test compared to all tissue-specific genes:")
result<-wilcox.test(ts_dnds$dN.dS, ov_dnds$dN.dS)
print(result)

#Li et al testis specific
testis_prots0<-ligenes$protID[which(ligenes$tissue=="Te")]
testis_prots<-unlist(sapply(testis_prots0, function(x) unlist(strsplit(x, split=";"))))

print(head(testis_prots))
print(length(testis_prots))

testis_dnds<-dnds[which(dnds$protID %in% testis_prots),]
print(paste("Number of tissue-specific proteins in testis:", nrow(testis_dnds)))
print(paste("Median dN/dS for testis-specific proteins:", median(testis_dnds$dN.dS, na.rm=TRUE)))
print("Wilcox test compared to all genes:")
result<-wilcox.test(dnds$dN.dS, testis_dnds$dN.dS)
print(result)
print("Wilcox test compared to all tissue-specific genes:")
result<-wilcox.test(ts_dnds$dN.dS, testis_dnds$dN.dS)
print(result)
print("Wilcox test compared to ovary-specific genes:")
result<-wilcox.test(ov_dnds$dN.dS, testis_dnds$dN.dS)
print(result)

#Chalmel et al. testis-specific
chalmel_ts<-scan("/projects/ek112884e/murines_revisions/FROM_REDWOOD/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_testisSpecific.txt", what=character())
chalmel_prots<-intersect(chalmel_ts, prots)

print(head(chalmel_prots))
print(length(chalmel_prots))
print(paste("Overlap between Li et al. and Chalmel et al. testis-specific lists:", length(intersect(chalmel_prots, testis_prots))))

chalmel_dnds<-dnds[which(dnds$protID %in% chalmel_prots),]
print(paste("Number of tissue-specific proteins in chalmel:", nrow(chalmel_dnds)))
print(paste("Median dN/dS for chalmel-specific proteins:", median(chalmel_dnds$dN.dS, na.rm=TRUE)))

print(paste("Overlap between Li et al. and Chalmel et al. testis-specific lists (prots in our dataset only):", length(intersect(testis_dnds$protID, chalmel_dnds$protID))))

print("Wilcox test compared to all genes:")
result<-wilcox.test(dnds$dN.dS, chalmel_dnds$dN.dS)
print(result)
print("Wilcox test compared to all tissue-specific genes:")
result<-wilcox.test(ts_dnds$dN.dS, chalmel_dnds$dN.dS)
print(result)
print("Wilcox test compared to ovary-specific genes:")
result<-wilcox.test(ov_dnds$dN.dS, chalmel_dnds$dN.dS)
print(result)
