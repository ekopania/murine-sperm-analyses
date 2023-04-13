#PURPOSE: Run permulations on GO enrichment test
#Modified from getPermsBinary()
	#I copied the enrichment parts of this function and modified them to work with topGO, instead of MSigDB

#If I decide to use MSigDB annotations, here they are for mouse: https://bioinf.wehi.edu.au/MSigDB/v7.1/

#Calls functions in this script
source("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/customPermulationScripts.r")

#' @param  perms   Result from a permulations run containing perm p vals and perm stats
#' @param  GO_thresh  The threshold cutoff to test for GO enrichment (number of genes)
#' @return data  A list object with P values and stat values for each permulation's topGO analysis; stat values are observed / expected number of highly accelerated genes in the annotation set
getGOEnrichPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, method="k", min.pos=2, perms, GO_thresh){
    require(topGO)
    #Run GO analysis for real data
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat = sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    print(head(realstat))
    realenrich = accGO(realstat, GO_thresh)

    print(head(realenrich))

    #get permulation stat from permulation run; should correspond to corStat from permulation output
    #make sure loci are in the same order as the real data
    permStatvals = perms$corStat[match(names(realstat), rownames(perms$corStat)),]

    realenrich=realenrich[order(realenrich$GO.ID),] 
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich)))
    rownames(newdf)=realenrich$GO.ID
    permenrichP=newdf
    permenrichStat=newdf
    
    counter=1;
    while (counter <= numperms){
      print(paste("Working on perm", counter))
      stat = permStatvals[,counter]
      print(head(stat))
      names(stat) = rownames(permStatvals)
      print(head(stat))
      enrich=accGO(na.omit(stat), GO_thresh)
      #sort and store enrichment results
      current=enrich
      rownames(current)=current$GO.ID
      enrich=current[order(rownames(current)),]
      enrich=enrich[match(rownames(permenrichP), rownames(enrich)),]
      permenrichP[,counter]=enrich$pval
      permenrichStat[,counter]=enrich$stat
      counter=counter+1
    }

    data=vector("list", 2)
    data[[1]]=permenrichP
    data[[2]]=permenrichStat
    names(data)=c("enrichP", "enrichStat")

    data
}

#Helper function to convert between protein IDs and gene IDs using Ensembl biomaRt
#I wrote this because the murine data are in protID form but topGO annotations are in geneID form
#Uses the Nov 2020 version of Ensembl because that is what we have been using for murine stuff
protToGene=function(in_prots){
    require(biomaRt)
    ens_mus = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org/")
    all_genes = getBM(attributes=c('ensembl_peptide_id','ensembl_gene_id','external_gene_name'), mart=ens_mus)
    colnames(all_genes) = c("prot_id","gene_id","gene_name")
    out_genes<-c()
    for(i in in_prots){
      if(i %in% all_genes$prot_id){
        out_genes<-c(out_genes, all_genes$gene_id[which(all_genes$prot_id==i)])
      } else{
        out_genes<-c(out_genes, NA)
      }
    }

    out_genes 
}

#' Runs topGO based on an input RERconverge result and threshold cutoff for top most accelerated genes
#' @param  res  Stat from an RERconverge run (signed log10 p-value for each gene)
#' @param  thresh  Number of genes to include in the foreground set for topGO analysis (i.e., if you want to test for GO enrichment among the top 100 most accelerated genes, set thresh = 100)
#' @return  final_result  A table of topGO outputs, including information on GO terms, numbers of genes annotated, observed, and expected, p-values, and fold-change enrichment
accGO=function(res, thresh){
    require(topGO)
    #res_noNA<-res[which(!(is.na(names(res))))]
    #res_acc<-res_noNA[which(res_noNA > 0)]
    #print(sum(is.na(names(res_acc))))
    #Sort res such that the most accelerated genes will come first in the vector
    res_sorted<-res[order(res)]
    #Get top most accelerated genes
    subset_acc<-names(res_sorted[order(res_sorted, decreasing=TRUE)])[1:thresh]
    print(head(res_sorted[order(res_sorted, decreasing=TRUE)]))
    print(paste("Number of genes in focal accelerated subset:", length(subset_acc)))
    print(paste("Total number of genes in accelerated dataset:", length(res_sorted)))
    #print(head(subset_acc))
    total<-names(res_sorted)
    #Set up functions for topGO
    inSubset<-c()
    for(gene in total){
        if(gene %in% subset_acc){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
    }
    #print(length(subset_acc))
    #print(length(which(inSubset==1)))
    if(length(subset_acc) != length(which(inSubset==1))){
        stop("ERROR: length of 'subset_acc' file and 'inSubset' vector don't match up!")
    }
    names(inSubset)<-total
    geneSelFunc<-function(iO){
        return(iO==1)
    }
    #Run topGO and parse results
    myGOdata<-new("topGOdata", description="topAccelerated", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
    result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
    resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
    #Correct p-values for multiple tests
    resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
    #print(head(resultData))
    
    #Get fold enrichment stat
    foldEnrich<-resultData$Significant / resultData$Expected

    #Make final table to output (topGO output table, plus corrected p-value and fold-enrichment stat)
    final_result<-as.data.frame(cbind(resultData, stat=foldEnrich))
    colnames(final_result)<-c("GO.ID", "Term", "Num.Annotated.Genes", "Num.Significant.Genes", "Num.Expected.Genes", "pval", "p.adj", "stat")

    final_result
}

#' Get p-values for GO enrichment permulation results
#' @param  realenrich    The result of a topGO analysis with the real data from RERconverge
#' @param  permvals      The output of a permulation run
#' @param  progress      Boolean for whether or not to occassionally update on progress in the while loop
#' @return  enrichpvals  A list of permulation p-values
permpvalGOenrich=function(realenrich, permvals, progress=FALSE){
    current=realenrich
    rownames(current)=realenrich$GO.ID
    
    #Make sure real enrichment and permulation enrichment tables are in the same orders (sorted by GO ID)
    realenrich=current[match(rownames(permvals$enrichStat), rownames(current)),]

    permenrich=permvals$enrichStat
    enrichpvals=list()
    rowlen=nrow(permenrich)
    print(paste("Looping through", rowlen, "GO terms..."))
    rowcount=1
    pvallist=c()
    #Loop through each permulation
    while(rowcount<=rowlen){
        #Prints out progress to keep track of which permulation loop is on
        if( progress && ((rowcount %% 100)==0) ){
            print(paste("On row", rowcount))
        }
        #print(realenrich[rowcount,]$stat)
        #Checks if stat (fold-enrichment) is NA for this permulation; if so return NA for permulation p-val
        if(is.na(realenrich[rowcount,]$stat)){
          #pval=lessnum/denom
          pval=NA
        #Get the p-val for this permulation
        }else{
          #Counts the number of permulated fold-enrichment values that are greater than the real fold-enrichment value for this GO term
          lessnum=sum(abs(permenrich[rowcount,])>abs(realenrich[rowcount,]$stat), na.rm=T)
          #Counts the total number of non-NA permulation enrichment stats
          denom=sum(!is.na(permenrich[rowcount,]))
          #Add 1 to both numerator and denominator so that no perm p-vals are zero (this is how the normal permulation functions work)
          pval=(lessnum+1)/(denom+1)
        }
        pvallist=c(pvallist, pval)
        rowcount=rowcount+1
    }
    names(pvallist)=rownames(realenrich)
    enrichpvals=as.list(pvallist)
    
    enrichpvals
}
