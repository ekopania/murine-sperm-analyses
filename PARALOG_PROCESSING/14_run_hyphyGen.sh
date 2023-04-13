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
#Full paralog set
#BUSTED-PH
#python /home/ek112884/software/core/hyphy-interface/hyphy_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_MERGE_PARALOGS-f50-seq20-site50/nt/ -m busted-ph -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_MERGED_PARALOGS-f50-seq20-site50/loci/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_busted-ph-OUmodelBackground -tb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -p ~/software/hyphy-analyses/ -part good_lab_cpu -tasks 9 -cpus 1 -mem 0 -n busted-ph_top10_OUbackground
#SLAC
#python /home/ek112884/software/core/hyphy-interface/hyphy_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_MERGE_PARALOGS-f50-seq20-site50/nt/ -m slac -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_MERGED_PARALOGS-f50-seq20-site50/loci/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_slac -part good_lab_reincarnation -tasks 10 -cpus 1 -mem 0 -n slac_full_paralog_set

#RTM dataset
#BUSTED-PH
#python /home/ek112884/software/core/hyphy-interface/hyphy_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m busted-ph -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_busted-ph-OUmodelBackground -tb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -p ~/software/hyphy-analyses/ -part good_lab_cpu -tasks 9 -cpus 1 -mem 0 -n busted-ph_top10_OUbackground_RTM
#SLAC
#python /home/ek112884/software/core/hyphy-interface/hyphy_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m slac -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_slac  -part good_lab_reincarnation -tasks 10 -cpus 1 -mem 0 -n slac_RTM

#RTM dataset - BY EXON
#BUSTED-PH
python ~/software/core/hyphy-interface/hyphy_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m busted-ph -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci/ -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_busted-ph-OUmodelBackground -tb /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -part clarkn-shared-rw -tasks 20 -cpus 1 -mem 2000 -n busted-ph_top10_OUbackground_RTM_byExon
