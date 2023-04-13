#!/bin/bash
#SBATCH --job-name=paml_gen
#SBATCH --output=paml_gen-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1-0 #Set time limit to 7 days zero hours (overrides server default of 2 day limit)
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:

###PARALOG DATASET###
#M0
#python ~/software/core/paml-interface/paml_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_MERGE_PARALOGS-f50-seq20-site50/nt/ -m m0  -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_MERGED_PARALOGS-f50-seq20-site50/loci -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_m0 -part good_lab_reincarnation -tasks 9 -cpus 1 -mem 0 -n repro_paralogs_merged_paralogs_m0

#M1a vs M2a
#python ~/software/core/paml-interface/paml_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_MERGE_PARALOGS-f50-seq20-site50/nt/ -m m1a  -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_MERGED_PARALOGS-f50-seq20-site50/loci -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_m1a -part good_lab_reincarnation -tasks 9 -cpus 1 -mem 0 -n repro_paralogs_merged_paralogs_m1a
#python ~/software/core/paml-interface/paml_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED_MERGE_PARALOGS-f50-seq20-site50/nt/ -m m2a  -genetrees /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE_MERGED_PARALOGS-f50-seq20-site50/loci -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SELEC_TESTS/MERGED_PARALOGS-f50-seq20-site50_m2a -part good_lab_reincarnation -tasks 9 -cpus 1 -mem 0 -n repro_paralogs_merged_paralogs_m2a

###PARALOG DATASET - RTM ONLY###
##M0
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m m0 -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_m0 -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_m0
#
##M1a vs M2a
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m m1a -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_m1a -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_m1a
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m m2a -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_m2a -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_m2a
#
##BS test
##Branch-site test (BS2)
##See PAML manual pg 30-31
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m bs -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_branchSite -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_branchSite
##Null model for comparison for branch-site test (BS1)
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt/ -m bs_null -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50/loci -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM-f16-seq20-site50_branchSiteNull -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_branchSiteNull

###PARALOG DATASET - RTM ONLY - BY EXON - PARALOGS REMOVED###
#M0
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m m0 -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m0 -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_byExon_m0
#Remove second part of sample names from fasta and tree, which PAML doesn't like
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m0/*filter/codeml.fa | while read file; do sed -i 's/_|_.*//' ${file}; done
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m0/*filter/codeml.tre | while read file; do sed -i 's/_|_Mus_musculus_[[:alnum:]]\+-p[[:digit:]]\+:/:/g' ${file}; done


#M1a vs M2a
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m m1a -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m1a -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_byExon_m1a
#python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m m2a -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m2a -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_byExon_m2a
#Remove second part of sample names from fasta and tree, which PAML doesn't like
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m1a/*filter/codeml.fa | while read file; do sed -i 's/_|_.*//' ${file}; done
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m1a/*filter/codeml.tre | while read file; do sed -i 's/_|_Mus_musculus_[[:alnum:]]\+-p[[:digit:]]\+:/:/g' ${file}; done
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m2a/*filter/codeml.fa | while read file; do sed -i 's/_|_.*//' ${file}; done
#ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m2a/*filter/codeml.tre | while read file; do sed -i 's/_|_Mus_musculus_[[:alnum:]]\+-p[[:digit:]]\+:/:/g' ${file}; done

#BS test
#Branch-site test (BS2)
#See PAML manual pg 30-31
python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m bs -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_branchSite -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_byExon_branchSite
#Null model for comparison for branch-site test (BS1)
python ~/software/core/paml-interface/paml_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt-removed-paralogs/ -m bs_null -genetrees /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_byExon_RTM-f16-seq20-site50/loci -targetclade /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/large_testes_species.q90.txt -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_branchSiteNull -part clarkn-shared-rw -tasks 12 -cpus 1 -mem 0 -n repro_paralogs_RTM_byExon_branchSiteNull
