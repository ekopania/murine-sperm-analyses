#PURPOSE: Test if genes with RERs associated with trait are enriched for GO terms/bio functions

library(RERconverge)
library(biomaRt)

print("Running enrichment test for RERconverge output on binary trait")

print("Loading data and getting stats from RERconverge run...")
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")
stats<-getStat(res)

print("Reading annotation file and converting to correct format...")
annots<-read.gmt("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/c2.all.v7.4.symbols.gmt")
annotlist=list(annots)
names(annotlist)="MSigDBpathways"

print("Converting mouse protein names to human gene names...")
#USE unlist(sapply) and biomart to go from stat names (e.g., ENSMUSP00000000010-mafft-cds.filter.AAreplace) to human gene names (e.g., HOXB9)
mouse_prot<-unlist(sapply(names(stats), function(x) gsub("-.*","",x)))
#Using an old version of Ensembl because getLDS() has issues after the April 2022 update
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org/")
ens_hum<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "https://nov2020.archive.ensembl.org/")
all_genes<-getBM(attributes=c('ensembl_peptide_id','external_gene_name'), mart=ens_mus)
colnames(all_genes)<-c("prot_id","gene_name")
musNames<-unlist(sapply(mouse_prot, function(x) all_genes$gene_name[which(all_genes$prot_id==x)]))
hum_genes<-getLDS(attributes=c("ensembl_peptide_id","external_gene_name"), values=musNames, mart=ens_mus, attributesL=c("external_gene_name"), martL=ens_hum, uniqueRows=T)
colnames(hum_genes)<-c("mus_ID", "mus_name", "hum_name")
print(dim(hum_genes))
print(head(hum_genes))
names(stats)<-mouse_prot
print(head(stats))

#Some stuff to deal with 1:many and many:1 orthos, as well as genes that don't have an ortholog in both species
#If many human genes, keep the first one
#If many mouse genes, keep the last one
#If no ortholog, remove from dataset
hum_genes_sorted<-hum_genes[order(hum_genes$mus_ID),]
hum_genes_uniq<-hum_genes_sorted[which(!(duplicated(hum_genes_sorted$mus_ID))),]
hum_genes_final<-hum_genes_uniq[which(hum_genes_uniq$mus_ID != ""),]

stats_humNames<-c()
no_ortho<-c()
multi_mouse<-c()
for(i in 1:length(stats)){
	if(names(stats)[i] %in% hum_genes_final$mus_ID){
		hum_g<-hum_genes_final$hum_name[which(hum_genes_final$mus_ID==names(stats)[i])]
			if(hum_g %in% names(stats_humNames)){
				#print(paste("Two mouse genes for the following human gene:", hum_g))
				multi_mouse<-c(multi_mouse, names(stats)[i])
				stats_humNames<-c(stats_humNames, stats[i])
				names(stats_humNames)[length(stats_humNames)]<-hum_g
			} else{
				stats_humNames<-c(stats_humNames, stats[i])
				names(stats_humNames)[length(stats_humNames)]<-hum_g
			}
	} else{
		no_ortho<-c(no_ortho, names(stats)[i])
	}
}
print(paste("Number of genes removed from stats due to no ortholog:", length(no_ortho)))
print(paste("Number of genes with multiple mouse genes corresponding to a single human ortholog:", length(multi_mouse)))
print(length(stats_humNames))
print(head(stats_humNames))

print("Running enrichment test...")
enrichment<-fastwilcoxGMTall(stats_humNames, annotlist, outputGeneVals=T, num.g=10)
enrichment_sorted<-enrichment$MSigDBpathways[order(enrichment$MSigDBpathways$p.adj),]

print(head(enrichment_sorted[,c("stat","pval","p.adj","num.genes")]))

#All the gene sets I could find that seemed related to reproduction
repro_sets<-c("MATZUK_SPERMATID_DIFFERENTIATION", "MATZUK_SPERMATOCYTE", "MATZUK_SPERMATOGONIA", "MATZUK_SPERMATOZOA", "REACTOME_SPERM_MOTILITY_AND_TAXES", "WEBER_METHYLATED_HCP_IN_SPERM_DN", "WEBER_METHYLATED_HCP_IN_SPERM_UP", "WEBER_METHYLATED_ICP_IN_SPERM_DN", "WEBER_METHYLATED_ICP_IN_SPERM_UP", "WEBER_METHYLATED_LCP_IN_SPERM_DN", "WEBER_METHYLATED_LCP_IN_SPERM_UP", "CHEN_ETV5_TARGETS_TESTIS", "CHEN_ETV5_TARGETS_SERTOLI", "REACTOME_TRANSCRIPTIONAL_REGULATION
+ _OF_TESTIS_DIFFERENTIATION", "SU_TESTIS", "YOKOE_CANCER_TESTIS_ANTIGENS", "MATZUK_MALE_REPRODUCTION_SERTOLI", "MATZUK_CENTRAL_FOR_FEMALE_FERTILITY", "MATZUK_CUMULUS_EXPANSION", "MATZUK_EARLY_ANTRAL_FOLLICLE", "MATZUK_FERTILIZATION", "MATZUK_IMPLANTATION_AND_UTERINE", "MATZUK_LUTEAL_GENES", "MATZUK_MEIOTIC_AND_DNA_REPAIR", "MATZUK_OVULATION", "MATZUK_POSTIMPLANTATION_AND_POSTPAR
+ TUM", "MATZUK_PREOVULATORY_FOLLICLE","MATZUK_STEROIDOGENESIS","REACTOME_REPRODUCTION", "WP_MALE_INFERTILITY")
print(paste("Number of repro gene sets in enrichment result:", length(which(repro_sets %in% rownames(enrichment_sorted)))))
print("Enrichment results for repro sets:")
print(enrichment_sorted[which(rownames(enrichment_sorted) %in% repro_sets),c("stat","pval","p.adj","num.genes")])

#save data
write.table(enrichment_sorted, file="RERconverge_output.logRTM_binary.enrichment.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
save(stats, stats_humNames, enrichment, enrichment_sorted, file="RERconverge_output.logRTM_binary.enrichment.rds")

print("Done with 11.2_enrichment.binary.r")
