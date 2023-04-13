#!/bin/bash
#PURPOSE: Use slurm to run Rscripts for l1ou
#
# Job name:
#SBATCH --job-name=l1ou
#SBATCH --output=l1ou-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@genetics.utah.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:
module load libxml2
module load R/4.1.1

#echo "Pruning..."
#Rscript 01_prune_tree.R

echo "Running l1ou..."
Rscript 02_run_l1ou.r

echo "Done!"
