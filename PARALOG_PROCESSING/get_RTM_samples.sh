#!/bin/bash
#PURPOSE: Get gene fastas with only the species that have RTM data
#
# Job name:
#SBATCH --job-name=get_RTM_gene_fasta
#SBATCH --output=get_RTM_gene_fasta-%j.log
##SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
###SBATCH --mail-user=emily.kopania4@genetics.utah.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=4G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
#SBATCH --time=7-0 #Set time limit to 7 days zero hours (overrides server default of 2 day limit)
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
#
## Command(s) to run:
module load seqtk

ls /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/GENE_FASTAS/*fa | while read file; do
	name=$(echo "${file}" | cut -d "/" -f 8)
	echo "Working on ${name}"
	rm temp_samples.txt
	cat rtm_spec.txt | while read line; do #Need this to get the full sample names with names, paralog numbers, etc
		grep ">${line}" "${file}" >> temp_samples.txt
	done
	sed -i 's/>//' temp_samples.txt
	seqtk subseq "${file}" temp_samples.txt > "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/GENE_FASTAS_RTM/${name}"
done

#cleanup
rm temp_samples.txt

echo "Done!"
