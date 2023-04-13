#!/bin/bash
#PURPOSE: Get exon IDs that correspond to exome capture targeted regions using mouse GTF file
#
# Job name:
#SBATCH --job-name=get_targeted_exons
#SBATCH --output=get_targeted_exons-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=2G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#Use bedtools intersect to get all lines from the mouse GTF that overlap with targeted regions
#-wb prints the gtf lines as well as the bedfile lines
bedtools intersect -a mm9_to_mm10_targeted_regions.bed -b /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf -wb > mm10_targeted_regions.gtf

#Isolate gtf lines that have peptide IDs
grep "ENSMUSE" mm10_targeted_regions.gtf > mm10_targeted_regions.exonOnly.gtf

#Use sed to remove everything before the peptide ID starts
sed 's/^.*exon_id "//' mm10_targeted_regions.exonOnly.gtf > sed1.gtf

#Use sed to remove everything after the peptide ID ends
sed 's/"; .*//' sed1.gtf > exon_list.all_targets.txt
#Left with a list of exon IDs that overlap targeted regions (mm10)

#Clean up
rm mm10_targeted_regions.gtf
rm mm10_targeted_regions.exonOnly.gtf
rm sed1.gtf

echo "Done!" 
