#!/bin/bash
#SBATCH --job-name=hyphy_gen
#SBATCH --output=hyphy_gen-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

#Include hyphy scripts in path
export PATH=$PATH:/ix3/nclark/ekopania/MURINAE_REVISIONS/core/hyphy-interface/
echo $PATH

#Run hyphy_gen.py
#Use full tree for repro dataset
#Set up with synonymous rate variation and multi-nucleotide mutations
python /ix3/nclark/ekopania/MURINAE_REVISIONS/core/hyphy-interface/hyphy_gen.py -i /ix3/nclark/ekopania/MURINAE_REVISIONS/alns/ -m relax -genetrees /ix3/nclark/ekopania/MURINAE_REVISIONS/GENE_TREES/ -o /ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_relax/OUfg_fullReproSet_GOsperm -tb OU_fgs.txt -p hyphy -part good_lab_cpu -tasks 64 -cpus 1 -mem 12000 -n OUfg_fullReproSet_GOsperm #-rb OU_bgs.txt

#Use pared tree to speed up runtime (used this one for PAML analyses); use full coding set alignments
#python /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE_REVISIONS/core/hyphy-interface/hyphy_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE_REVISIONS/alns/ -m relax -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE_REVISIONS/GENE_TREES/ -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE_REVISIONS/HYPHY_relax/OUfg_fromFullCoding_noOutgroups -tb OU_fgs.txt -p hyphy -part clarkn-shared-rw -tasks 60 -cpus 1 -mem 0 -n relax_OUfg_fromFullCoding # -rb OU_bgs.txt

#TEST
#python TEMP_CORE/core/hyphy-interface/hyphy_gen.py -i TEST_ALN -m busted-ph -genetrees TEST_TREES -o /mnt/beegfs/ek112884/murinae/DNDS_BY_SPERM_STAGE/HYPHY/BUSTED-PH/TEST_OUTPUT -tb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -rb /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/small_testes_species.OUmodel.txt -p ~/software/hyphy-analyses/ -part good_lab_cpu -tasks 60 -cpus 1 -mem 0 -n busted-ph_RTM_TEST
