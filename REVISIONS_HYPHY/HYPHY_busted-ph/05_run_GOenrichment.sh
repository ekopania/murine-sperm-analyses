#!/bin/bash
#PURPOSE: Run 06_GOenrichment.r
#SBATCH --job-name=GOenrichment
#SBATCH --output=GOenrichment-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

source activate /ihome/nclark/emk270/software/envs/r4.1.1

Rscript 06_GOenrichment.r
