#PURPOSE: Make a giant table combining RERconverge, dNdS, and selection test results

#Read in molecular evolution data; format protein ID names if necessary
rerconverge<-read.csv("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/RERconverge_output.logRTM_binary_residuals.OUmodel.RTMspeciesOnly.withPermulationsSSM.customSSM.customSisters.csv", header=TRUE)
rerconverge$X<-unlist(sapply(rerconverge$X, function(x) sub("-mafft-cds.filter.AAreplace", "", x)))

dnds<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/RTM_SET_pared_M0/dnds.RTM_SET_pared_M0.txt", header=TRUE)

posSelec<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.m1aVm2a.txt", header=TRUE)
#relax<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.m1aVbsnull.txt", header=TRUE)
relax<-read.table("/scratch/general/pe-nfs1/kopania/murines/PAML/pos_selec.pValCor.sigOnly.posOnly.m1aVbsnull.txt", header=TRUE)
bs<-read.table("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.pValCor.bsTest.txt", header=TRUE)

#Read in gene sets for tissues/cell types
bd<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BD.any.txt", what=character())
bu<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.BU.any.txt", what=character())
cg<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.CG.any.txt", what=character())
dp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.DP.any.txt", what=character())
vp<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.VP.any.txt", what=character())
sv<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.SV.any.txt", what=character())

som<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.somatic.txt", what=character())
spg<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatogonia.txt", what=character())
pre<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.prelep.txt", what=character())
spc<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatocytes.txt", what=character())
spd<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.spermatids.txt", what=character())
elo<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt", what=character())

ts<-scan("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/prot_list_testisSpecific.txt", what=character())


#Get set of protein-coding loci across all these
loci<-Reduce(union, list(rerconverge$X, dnds$protID, posSelec$protID, bs$protID, relax$protID))

#Get list of tissue and cell types sets protein is in
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

full_df<-c()
for(i in loci){
	temp_vec<-c(i, membership_lists[i])
	if(i %in% rerconverge$X){
		idx<-which(rerconverge$X==i)
		temp_vec<-c(temp_vec, rerconverge$Rho[idx], rerconverge$N[idx], rerconverge$P[idx], rerconverge$p.adj[idx], rerconverge$permpval[idx], rerconverge$permpvaladj[idx])
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA, NA)
	}
	if(i %in% dnds$protID){
		idx<-which(dnds$protID==i)
		temp_vec<-c(temp_vec, dnds$dN.dS[idx])
	} else{
		temp_vec<-c(temp_vec, NA)
	}
	if(i %in% posSelec$protID){
		idx<-which(posSelec$protID==i)
		if(posSelec$lnL[idx] >= posSelec$lnL.null[idx]){ #Make sure model converged
			temp_vec<-c(temp_vec, posSelec$lnL[idx], posSelec$lnL.null[idx], posSelec$diff[idx], posSelec$P[idx], posSelec$p.adj[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
	}
	if(i %in% bs$protID){
		idx<-which(bs$protID==i)
		if(bs$lnL[idx] >= bs$lnL.null[idx]){ #Make sure model converged
			temp_vec<-c(temp_vec, bs$lnL[idx], bs$lnL.null[idx], bs$diff[idx], bs$P[idx], bs$p.adj[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
	}
	if(i %in% relax$protID){
		idx<-which(relax$protID==i)
		if(relax$lnL[idx] >= relax$lnL.null[idx]){ #Make sure model converged
			temp_vec<-c(temp_vec, relax$lnL[idx], relax$lnL.null[idx], relax$diff[idx], relax$P[idx], relax$p.adj[idx])
		} else{
			temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
		}
	} else{
		temp_vec<-c(temp_vec, NA, NA, NA, NA, NA)
	}

	full_df<-rbind(full_df, temp_vec)
}

full_df<-as.data.frame(full_df)
print(dim(full_df))
colnames(full_df)<-c("protID", "repro_gene_sets", "RERconverge_Rho", "RERconverge_Nspecies", "RERconverge_P", "RERconverge_p.adj", "RERconverge_permpval", "RERconverge_permpval.adj", "dN.dS", "lnL.m2", "lnL.m1", "m1Vm2.lrt", "m1Vm2.P", "m1Vm2.p.adj", "bs.lnL", "bs.lnL.null", "bs.lrt", "bs.P", "bs.p.adj", "relax.lnL", "relax.lnL.null", "relax.lrt", "relax.P", "relax.p.adj")

write.csv(full_df, "all_molec_evo_results.csv", row.names=FALSE)
