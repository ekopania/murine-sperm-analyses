#!/bin/bash
##PURPOSE: Run Rscript 08_find_filtered_step.r from slurm
#
# Job name:
#SBATCH --job-name=find_filtered
#SBATCH --output=find_filtered-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=20G #Not sure if I should mess with these...
#SBATCH --mem=24G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate r4

ls prot_list.*any*txt prot_list.*strict*txt prot_list.greenEtal2018.*txt | while read file; do
	echo "${file}"
	myname=$(echo "${file}" | cut -d "." -f 2-3)
	echo "Working on ${myname}"
	Rscript 08_find_filtered_step.r "${file}" > "filter_stats.${myname}.txt"
done

echo "Done!"
