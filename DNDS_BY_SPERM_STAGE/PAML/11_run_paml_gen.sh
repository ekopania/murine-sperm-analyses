#!/bin/bash
#SBATCH --job-name=paml_gen
#SBATCH --output=paml_gen-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100
#SBATCH --time=1-0 #Set time limit to 7 days zero hours (overrides server default of 2 day limit)
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:

##BS test
##Branch-site test (BS2)
##See PAML manual pg 30-31
python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/repro_gene_alignments/ -m bs -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/repro_trees -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o RTM_SET_pared_branchSite_reproGenes_top10 -part clarkn-shared-rw -tasks 10 -cpus 1 -mem 2000 -n reproGenes_RTM_pared_branchSite_top10
##Null model for comparison for branch-site test (BS1)
python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/repro_gene_alignments/ -m bs_null -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/repro_trees -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o RTM_SET_pared_branchSiteNull_reproGenes_top10 -part clarkn-shared-rw -tasks 10 -cpus 1 -mem 2000 -n reproGenes_RTM_pared_branchSiteNull_top10
