#!/bin/bash
#SBATCH --job-name=hyphy_parse
#SBATCH --output=hyphy_parse-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

#Include hyphy scripts in path
export PATH=$PATH:/ihome/nclark/emk270/software/core/hyphy-interface/
#echo $PATH

#Run hyphy_to_csv.py - SRV and MNM 2+3
#python /ix3/nclark/ekopania/MURINAE_REVISIONS/core/hyphy-interface/hyphy_to_csv.py -i OUfg_fullReproSet/ -m busted -o busted_output.OUfg_fullReproSet.csv

#Run hyphy_to_csv.py - no SRV or MNMs
python /ihome/nclark/emk270/software/core/hyphy-interface/hyphy_to_csv.py -i OUfg_fullReproSet_noSRV_noMNM/ -m busted -o busted_output.OUfg_fullReproSet.noSRV_noMNM.csv -mns None

#Run hyphy_to_csv.py - yes SRV but no MNMs
#python /ihome/nclark/emk270/software/core/hyphy-interface/hyphy_to_csv.py -i OUfg_fullReproSet_yesSRV_noMNM/ -m busted -o busted_output.OUfg_fullReproSet.yesSRV_noMNM.csv -mns None

#Run hyphy_gen.py
#Use full tree for repro dataset
#python /home/ek112884e/software/core/hyphy-interface/hyphy_gen.py -i /mnt/beegfs/gt156213e/murinae-seq/03-Alignments/03-Exonerate/reproductive-testes-mass-coding-trimmed-f0-seq20-site50/nt/ -m relax -genetrees /mnt/lepus/home/gt156213e/murinae-seq/04-Phylo/reproductive-datasets/reproductive-testes-mass-coding-iqtree/loci -o /mnt/rattus/ek112884e/murinae/DNDS_BY_SPERM_STAGE/HYPHY/RELAX/OUfg_fullReproSet -tb OU_fgs.txt -rb OU_bgs.txt -p ~/software/hyphy-analyses/ -part good_lab_cpu -tasks 60 -cpus 1 -mem 0 -n relax_OUfg_fullReproSet

#Use pared tree to speed up runtime (used this one for PAML analyses); use full coding set alignments
#python3 /hellgate/home/ek112884e/software/core/hyphy-interface/hyphy_gen.py -i /projects/gt156213e/murinae-seq/03-Alignments/03-Exonerate/full-coding-trimmed-f175-seq20-site50/nt/ -m relax -genetrees /projects/ek112884e/murines_revisions/TREE_PARING_OLD/pared_tree.pruned_ultrametric_RTM.keepOUshifts.bl100.tree/GENE_TREES/ -o /projects/ek112884e/murines_revisions/HYPHY_RELAX/OUfg_paredSet_noOutgroups -tb OU_fgs.txt -p /hellgate/home/ek112884e/software/hyphy/hyphy -part good_lab_cpu -tasks 60 -cpus 1 -mem 0 -n relax_OUfg_paredSet # -rb OU_bgs.txt

#TEST
#python TEMP_CORE/core/hyphy-interface/hyphy_gen.py -i TEST_ALN -m busted-ph -genetrees TEST_TREES -o /mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/TEST_OUTPUT -tb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -p ~/software/hyphy-analyses/ -part good_lab_cpu -tasks 60 -cpus 1 -mem 0 -n busted-ph_RTM_TEST
