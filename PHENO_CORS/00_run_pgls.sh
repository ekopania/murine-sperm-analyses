#!/bin/bash
#PURPOSE: Use slurm to run Rscripts for pgls analysis
#
# Job name:
#SBATCH --job-name=pgls
#SBATCH --output=pgls-%j.log
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=4G #Not sure if I should mess with these...
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:
module load R/4.1.1

Rscript 01_pgls.r

echo "Done!"
