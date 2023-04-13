#!/bin/bash
#PURPOSE: Use bedtools to extract exon fasta sequences from mouse reference
#
# Job name:
#SBATCH --job-name=get_exon_fasta
#SBATCH --output=get_exon_fasta-%j.log
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

#echo "Formatting bed"
#sed -i 's/^/chr/' "../repro_exons.bed"
#sed -i 's/	-1	/	-	/' "../repro_exons.bed"
#sed -i 's/	1	/	+	/' "../repro_exons.bed"

#echo "Running bedtools getfasta"
#bedtools getfasta -name -fo "../repro_exons.fa" -fi "/mnt/beegfs/ek112884/REFERENCE_DIR/mm10.fa" -bed "../repro_exons.bed"

echo "Splitting fasta by exon"
faidx --split-files "../repro_exons.fa"
mv ENSMUSE*fa ../EXON_FASTAS


####OLD - Don't need this anymore because using Gregg's exon trimming script to get mm10 exons trimmed and in frame
#echo "Translating"
##cat ../repro_exons.bed | while read line; do
#	chr=$(echo "${line}" | cut -f 1)
#	sta=$(echo "${line}" | cut -f 2)
#	sto=$(echo "${line}" | cut -f 3)
#	exon=$(echo "${line}" | cut -f 4)
#	phase=$(echo "${line}" | cut -f 6)
#	phase1=$((${phase} + 1)) #AMAS starts at 1 but Ensembl starts at 0 for first nucleotide
#	echo "${chr} ${sta} ${sto} ${exon} ${phase1}"
#	if [ ${phase1} -eq 0 ]; then #-1 (0 after adding 1 for AMAS) are untranslated regions; skip
#		echo "Untranslated; skipping..."
#	else
#		python ~/software/AMAS/amas/AMAS.py translate -f fasta -d dna -i "../EXON_FASTAS/${exon}${chr}${sta}-${sto}.fa" --reading-frame ${phase1} --out-format fasta
#		mv "translated_${exon}${chr}${sta}-${sto}.fa-out.fas-out.fas" "../EXON_FASTAS"
#	fi
#done

echo "Done!"
