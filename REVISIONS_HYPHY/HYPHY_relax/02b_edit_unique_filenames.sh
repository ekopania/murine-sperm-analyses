#!/bin/bash
#SBATCH --job-name=unique_fnames
#SBATCH --output=unique_fnames-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

rm jobs/OUfg_fullReproSet_CPU1.uniques.sh
cat jobs/OUfg_fullReproSet_CPU1.sh | while read line; do
	#Copy over comment lines
	if [ $(echo ${line} | cut -c 1) == "#" ]; then
                echo ${line} >> jobs/OUfg_fullReproSet_CPU1.uniques.sh
	else
		loci=$(echo ${line} | cut -d "/" -f 8 | cut -d "-" -f 1)
		#If there were duplicate sequences, change the command so we only perform relaxh on the unique sequences and corresponding tree
        	if [ -f OUfg_fullReproSet_CPU1/${loci}-mafft-cds.filter/uniques.tre ]; then
                	newline=$(echo ${line} | sed -r 's/hyphyCPU/hyphy CPU/' | sed -r 's/relax.fa/uniques_noTree.fa/g' | sed -r 's/relax.tre/uniques.tre/g')
	        else
        	        newline=$(echo ${line} | sed -r 's/hyphyCPU/hyphy CPU/')
	        fi
        	echo ${newline} >> jobs/OUfg_fullReproSet_CPU1.uniques.sh
	fi
done
