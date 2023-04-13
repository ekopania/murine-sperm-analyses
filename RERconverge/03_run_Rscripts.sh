#!/bin/bash
#PURPOSE: Use slurm to run Rscripts for RERconverge
#
# Job name:
#SBATCH --job-name=RERconverge
#SBATCH --output=RERconverge-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@genetics.utah.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem=4G #Not sure if I should mess with these...
#SBATCH --time=7-0 #Set time limit to 7 days zero hours (overrides server default of 2 day limit)
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:
module load R/4.1.1

#echo "Pruning..."
#Rscript 04_prune_tree.R

#echo "Generating trees..."
#Rscript 05_make_trees.R

#echo "Running RERconverge"
#Rscript 07.1_get_RERs_logPheno.continuous.r #Log-transformed pheno data; continuous
#Rscript 07.2_get_RERs_logPheno.binary.r #Log-transformed pheno data; RTM as a binary trait

#echo "Running permulations..."
#Rscript 09.1_permulations_only.continuous.r
#Rscript 09.2_permulations_only.binary.r
#Rscript 09.3_permulations_only.binary.SSM.r

echo "Merging permulations..."
Rscript 10_combine_perms.r

#echo "Running enrichment"
#Rscript 11.1_enrichment.continuous.r
#Rscript 11.2_enrichment.binary.r

#echo "Running permulations on enrichment results..."
#Rscript 12.2_permulationsWenrichment.binary.R

echo "Running GO enrichment test..."
Rscript 13_GOenrichment.r

#echo "Running enrichment perms..."
#Rscript 16_enrichmentPermulation.r

#echo "Getting enrichment perm p values..."
#Rscript 17_getEnrichmentPermPs.r

echo "Done!"
