#!/bin/bash
#SBATCH --job-name=hyphy_toCSV
#SBATCH --output=hyphy_toCSV-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --partition=good_lab_reincarnation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0

#Include hyphy scripts in path
export PATH=$PATH:/home/ek112884/software/core/hyphy-interface/
echo $PATH

#Run hyphy_to_csv.py
#FULL PARALOG SET
python /home/ek112884/software/core/hyphy-interface/hyphy_to_csv.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_slac -m slac -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_slac/full-paralog-slac.csv --overwrite

#RTM SET
python /home/ek112884/software/core/hyphy-interface/hyphy_to_csv.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_slac -m slac -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_slac/rtm-slac.csv --overwrite
