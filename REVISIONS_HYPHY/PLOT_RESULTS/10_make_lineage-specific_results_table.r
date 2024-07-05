#PURPOSE: Make a table reporting BUSTED-PH and RELAX results for spermatogenesis GO-term genes

relax<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_relax/relax_output.GOspermatogenesisProteins_OUfg_fullReproSet.csv", header=TRUE, comment.char="#")
ph<-read.csv("/ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted-ph/busted-ph_output.GOspermatogenesisProteins_OUfg_fullReproSet.csv", header=TRUE, comment.char="#")
ph$protID<-unlist(sapply(ph$file, function(x) gsub("-mafft-cds.filter.json", "", x)))

loci<-union(relax$id, ph$protID)

print("Getting molecular evolution results...")
full_df<-c()
for(i in loci){
        temp_vec<-c(i)
        if(i %in% relax$id){
                idx<-which(relax$id==i)
		if(relax$LRT[idx] >= 0){ #Make sure model converged
        	#Header: id,LRT,p-value,K 
		temp_vec<-c(temp_vec, relax$LRT[idx], relax$p.value[idx], relax$K[idx])
        } else{
                temp_vec<-c(temp_vec, NA, NA, NA)
        }
        if(i %in% ph$protID){
                idx<-which(ph$protID==i)
                #Header: file,test_lrt,test_pval,bg_lrt,bg_pval,comp_lrt,comp_pval
		temp_vec<-c(temp_vec, ph$test_lrt[idx], ph$test_pval[idx], ph$bg_lrt[idx], ph$bg_pval[idx], ph$comp_lrt[idx], ph$comp_pval[idx])
	}
	} else{
                temp_vec<-c(temp_vec, NA, NA, NA, NA, NA, NA)
        }
	full_df<-rbind(full_df, temp_vec)
}

print("Making final table...")
full_df<-as.data.frame(full_df)
print(dim(full_df))
colnames(full_df)<-c("protID", "relax_lrt", "relax_P", "relax_K", "busted-ph_test_lrt", "busted-ph_test_P", "busted-ph_bg_lrt", "busted-ph_bg_P", "busted-ph_comparison_lrt", "busted-ph_comparison_P")

write.csv(full_df, "/ix3/nclark/ekopania/MURINAE_REVISIONS/PLOT_MOLEC_EVO/GOspermatogenesis_molec_evo_results.branch_site.csv", row.names=FALSE)
