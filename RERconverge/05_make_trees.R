setwd("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/")
source("estimateTreeFuncs.R")
#All alignments at once
#my_alndir<-"alignments_proteinCoding_trimmed_replaceAAs"
#Subset into directories to parallelize
my_alndir<-"alignments_proteinCoding_trimmed_replaceAAs_PARALOGS"
#Tree pruned to only include RTM data
my_intree<-"NEW_RTM_DATA_NOV2022/pruned_tree.RTM_taxa.tree"
my_output<-"NEW_RTM_DATA_NOV2022/trees_for_RERconverge_PARALOGS"
Sys.time()
estimatePhangornTreeAll(alndir=my_alndir, treefile=my_intree, output.file=my_output)
Sys.time()
