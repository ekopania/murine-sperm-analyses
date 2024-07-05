#PURPOSE: Make a giant table combining RERconverge, dNdS, and selection test results

#Read in molecular evolution data; format protein ID names if necessary
rerconverge<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/RERconverge_noParalogs/RERconverge_output.logRTM_binary_residuals.OUmodel.RTMspeciesOnly.withPermulationsSSM.customSSM.fromMasterTree.csv", header=TRUE)
rerconverge$X<-unlist(sapply(rerconverge$X, function(x) sub("-mafft-cds.filter.AAreplace", "", x)))
#Remove paralogs from RERconverge results
paralogs<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/RERconverge_noParalogs/paralog_loci.txt", what=character())
rerconverge<-rerconverge[which(!(rerconverge$X %in% paralogs)),]

#Read in data from all BUSTED runs
srv_dtMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.yesSRV_dtMNM.csv", header=TRUE, comment.char="#")
srv_noMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.yesSRV_noMNM.csv", header=TRUE, comment.char="#")
noSRV_noMNM<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted/busted_output.OUfg_fullReproSet.noSRV_noMNM.csv", header=TRUE, comment.char="#")

srv_dtMNM$protID<-unlist(sapply(srv_dtMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
srv_noMNM$protID<-unlist(sapply(srv_noMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))
noSRV_noMNM$protID<-unlist(sapply(noSRV_noMNM$file, function(x) gsub("-mafft-cds.filter.json", "", x)))

#Read in BUSTED model-averaged p-vals
model_avg_ps<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/PLOT_MOLEC_EVO/busted_output.OUfg_fullReproSet.model_averaged_pvals.csv", header=TRUE, col.names=c("protID","pval"))
model_avg_ps$p.adj<-p.adjust(model_avg_ps$pval, method="BH")

#Read in gene sets for tissues/cell types
som<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

ts<-scan("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.testisSpecific.txt", what=character())

bd<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BD.any.txt", sep="."), what=character())
bu<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.BU.any.txt", sep="."), what=character())
cg<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.CG.any.txt", sep="."), what=character())
dp<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.DP.any.txt", sep="."), what=character())
vp<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.VP.any.txt", sep="."), what=character())
sv<-scan(paste("/ix3/nclark/ekopania/MURINAE_REVISIONS/REPRO_PROTEIN_LISTS/prot_list.SV.any.txt", sep="."), what=character())


#Get set of protein-coding loci across all these
loci<-Reduce(union, list(rerconverge$X, srv_dtMNM$protID, srv_noMNM$protID, noSRV_noMNM$protID))

#Get list of tissue and cell types sets protein is in
print("Getting tissue and cell types...")
membership_lists<-c()
for(i in loci){
	groups<-c()
	if(i %in% bd){
		groups<-c(groups, "bulbourethral diverticulum")
	}
	if(i %in% bu){
		groups<-c(groups, "bulbourethral gland")
	}
	if(i %in% cg){
                groups<-c(groups, "coagulating gland")
        }
	if(i %in% dp){
                groups<-c(groups, "dorsolateral prostate")
        }
	if(i %in% vp){
                groups<-c(groups, "ventral prostate")
        }
	if(i %in% sv){
                groups<-c(groups, "seminal vescicle")
        }
	if(i %in% ts){
                groups<-c(groups, "testis-specific")
        }
	if(i %in% som){
                groups<-c(groups, "somatic")
        }
	if(i %in% spg){
                groups<-c(groups, "spermatogonia")
        }
	if(i %in% pre){
                groups<-c(groups, "pre-leptotene")
        }
	if(i %in% spc){
                groups<-c(groups, "spermatocytes")
        }
	if(i %in% spd){
                groups<-c(groups, "spermatid")
        }
	if(i %in% elo){
                groups<-c(groups, "elongating")
        }
	if(length(groups) > 0){
		all_groups<-paste(groups, collapse=";")
	} else{
		all_groups<-"none"
	}
	membership_lists<-c(membership_lists, all_groups)
}
names(membership_lists)<-loci

print("Getting molecular evolution results...")
full_df<-c()
for(i in loci){
	temp_vec<-c(i, membership_lists[i])
	if(i %in% rerconverge$X){
		idx<-which(rerconverge$X==i)
		temp_vec<-c(temp_vec, rerconverge$Rho[idx], rerconverge$N[idx], rerconverge$P[idx], rerconverge$p.adj[idx], rerconverge$permpval[idx], rerconverge$permpvaladj[idx])
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA, NA)
	}
	if(i %in% srv_dtMNM$protID){
		idx<-which(srv_dtMNM$protID==i)
		if(srv_dtMNM$lrt[idx] >= 0){ #Make sure model converged
			#Header: file,mnm2,mnm3,dn/ds,lrt,pval,aic
			temp_vec<-c(temp_vec, srv_dtMNM$dn.ds[idx], srv_dtMNM$lrt[idx], srv_dtMNM$pval[idx], srv_dtMNM$aic[idx], srv_dtMNM$mnm2[idx], srv_dtMNM$mnm3[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA, NA)
	}
	if(i %in% srv_noMNM$protID){
		idx<-which(srv_noMNM$protID==i)
		if(srv_noMNM$lrt[idx] >= 0){ #Make sure model converged
			temp_vec<-c(temp_vec, srv_noMNM$dn.ds[idx], srv_noMNM$lrt[idx], srv_noMNM$pval[idx], srv_noMNM$aic[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA)
	}
	if(i %in% noSRV_noMNM$protID){
		idx<-which(noSRV_noMNM$protID==i)
		if(noSRV_noMNM$lrt[idx] >= 0){ #Make sure model converged
			temp_vec<-c(temp_vec, noSRV_noMNM$dn.ds[idx], noSRV_noMNM$lrt[idx], noSRV_noMNM$pval[idx], noSRV_noMNM$aic[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA)
	}
	if(i %in% model_avg_ps$protID){
		idx<-which(model_avg_ps$protID==i)
		temp_vec<-c(temp_vec, model_avg_ps$pval[idx], model_avg_ps$p.adj[idx])
	} else{
		temp_vec<-c(temp_vec, NA, NA)
	}

	full_df<-rbind(full_df, temp_vec)
}

print("Making final table...")
full_df<-as.data.frame(full_df)
print(dim(full_df))
colnames(full_df)<-c("protID", "repro_gene_sets", "RERconverge_Rho", "RERconverge_Nspecies", "RERconverge_P", "RERconverge_p.adj", "RERconverge_permpval", "RERconverge_permpval.adj", "SRV_dtMNS_dN.dS", "SRV_dtMNS_lrt", "SRV_dtMNS_pval", "SRV_dtMNS_aic", "SRV_dtMNS_mns2", "SRV_dtMNS_mns3", "SRV_noMNS_dN.dS", "SRV_noMNS_lrt", "SRV_noMNS_pval", "SRV_noMNS_aic", "noSRV_noMNS_dN.dS", "noSRV_noMNS_lrt", "noSRV_noMNS_pval", "noSRV_noMNS_aic", "model_avg_P", "model_avg_p.adj")

write.csv(full_df, "/ix3/nclark/ekopania/MURINAE_REVISIONS/PLOT_MOLEC_EVO/all_molec_evo_results.csv", row.names=FALSE)
