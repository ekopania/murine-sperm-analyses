#!/bin/bash
#SBATCH --job-name=slac_branchAvgs
#SBATCH --output=slac_branchAvgs-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --partition=good_lab_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=0

source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro

#Include hyphy scripts in path
export PATH=$PATH:/home/ek112884/software/core/hyphy-interface/
echo $PATH

#Run branch_avgs.py
#FULL PARALOG SET - all genes
#python /home/ek112884/software/core/hyphy-interface/branch_avgs.py -t astral -n 10 -i /mnt/beegfs/gt156213e/murinae-seq/05-MolEvol/full-coding-astral-cf-rooted-rates.csv -r /mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/FULL_CODING/csv/ -f /mnt/beegfs/gt156213e/murinae-seq/05-MolEvol/full-coding-mg94-local-ds-filter-0.95quant.csv -s prot_list_all.full_coding_set.txt -o /mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/HYPHY/SLAC/full-coding-astral-rooted-rates.csv

#SPERM IMAGE SET - all genes
python /home/ek112884/software/core/hyphy-interface/branch_avgs.py -t astral -n 10 -i /mnt/beegfs/ek112884/murinae/reproductive-rtm-astral-rooted.csv -r /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_slac/csv/ -f /mnt/beegfs/gt156213e/murinae-seq/05-MolEvol/full-coding-mg94-local-ds-filter-0.95quant.csv -s prot_list_all.rtm_set.txt -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_slac/rtm-astral-rooted-rates.allGenes.csv 
