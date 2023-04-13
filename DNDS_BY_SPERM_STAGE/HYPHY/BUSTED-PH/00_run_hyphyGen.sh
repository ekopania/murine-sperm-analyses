#!/bin/bash
#SBATCH --job-name=hyphy_gen
#SBATCH --output=hyphy_gen-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --time=1-0 #Set time limit to 7 days zero hours (overrides server default of 2 day limit)
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:

#Include hyphy scripts in path
export PATH=$PATH:~/software/core/hyphy-interface/
echo $PATH

#Run hyphy_gen.py

#Repro genes dataset
#BUSTED-PH
python ~/software/core/hyphy-interface/hyphy_gen.py -i repro_gene_alignments/ -m busted-ph -genetrees repro_trees/ -o reproGenes_RTM_pared_top10-OUmodelBackground -tb /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -part clarkn-shared-rw -tasks 20 -cpus 1 -mem 2000 -n busted-ph_top10_OUbackground_reproGenes_RTM_pared
