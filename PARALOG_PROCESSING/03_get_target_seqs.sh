#!/bin/bash
#PURPOSE: Split fastas and re-combine them by gene rather than species
#
# Job name:
#SBATCH --job-name=get_target_seqs
#SBATCH --output=get_target_seqs-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=4G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro

#echo "Splitting species fastas"
#ls /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/paralogFinder/*/*fa | while read file; do
#	faidx --split-files "${file}"
#	mv *.fa ../SPLIT_FASTAS
#done

echo "Merging gene fastas"
cat ../repro_genes.txt | while read gene; do
	echo "${gene}"
	rm ../GENE_FASTAS/${gene}.fa #Rm if it already exists, instead of appending and possibly repeating some sequences
	for f in ../SPLIT_FASTAS/*${gene}*; do
		cat "${f}" >> "../GENE_FASTAS/${gene}.fa" #Append!
	done
done

echo "Done!"
