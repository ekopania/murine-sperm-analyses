#!/bin/bash
#PURPOSE: Split fastas and re-combine them by gene rather than species
#
# Job name:
#SBATCH --job-name=parallel_gen_get_target_seqs
#SBATCH --output=parallel_gen_get_target_seqs-%j.log
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

echo "Splitting species fastas"
echo "#!/bin/bash" > /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SPLIT_FASTAS/split_fastas_parallel.sh
ls /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/paralogFinder/*/*fa | while read file; do
	name=$(echo "${file}" | cut -d "/" -f 8 | cut -d "." -f 1)
	echo "${name}"
	#if [ -f /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SPLIT_FASTAS/${name}*.fa ]; then
	#	echo "Split fasta already exists for ${name}; skipping..."
	#else
		echo "faidx --split-files "${file}"" >> /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/SPLIT_FASTAS/split_fastas_parallel.sh
	#fi
done

echo "Done!"
