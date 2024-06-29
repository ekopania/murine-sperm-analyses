#!/bin/bash
#PURPOSE: Use hyphy functions to remove duplicate sequences from fasta and relabel pruned tree; should speed up computation
#https://github.com/veg/hyphy-analyses/tree/master/remove-duplicates
#https://github.com/veg/hyphy-analyses/tree/master/LabelTrees
#SBATCH --job-name=hyphy_rmdup
#SBATCH --output=hyphy_rmdup-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

source activate hyphy

for dir in /ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_busted-ph/OUfg_fullReproSet/*filter; do
	#Remove duplicate sequences
	hyphy /ihome/nclark/emk270/software/hyphy-analyses/remove-duplicates/remove-duplicates.bf --msa ${dir}/busted-ph.fa --tree ${dir}/busted-ph.tre --output ${dir}/uniques.fa ENV="DATA_FILE_PRINT_FORMAT=9"
	#Check if any duplicate sequences were removed
	if [ -f ${dir}/uniques.fa ]; then
		#Re-label the pruned output tree with TEST foreground branches
		hyphy /ihome/nclark/emk270/software/hyphy-analyses/LabelTrees/label-tree.bf --tree ${dir}/uniques.fa --label TEST --list /ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_relax/OU_fgs.sepLines.txt --output ${dir}/uniques.tre
		#Copy just the sequences (no tree) to a fasta file
		grep -v "(" ${dir}/uniques.fa  > ${dir}/uniques_noTree.fa
		sed -i '/^$/d' ${dir}/uniques_noTree.fa
	fi
done
