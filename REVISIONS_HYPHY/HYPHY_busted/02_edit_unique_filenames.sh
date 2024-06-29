#!/bin/bash
#SBATCH --job-name=unique_fnames
#SBATCH --output=unique_fnames-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

rm jobs/OUfg_fullReproSet_yesSRV_noMNM.uniques.sh
cat jobs/OUfg_fullReproSet_yesSRV_noMNM.sh | while read line; do
	#Copy over comment lines
	if [ $(echo ${line} | cut -c 1) == "#" ]; then
                echo ${line} >> jobs/OUfg_fullReproSet_yesSRV_noMNM.uniques.sh
	else
		loci=$(echo ${line} | cut -d "/" -f 8 | cut -d "-" -f 1)
		#If there were duplicate sequences, change the command so we only perform bustedh on the unique sequences and corresponding tree
        	if [ -f OUfg_fullReproSet_yesSRV_noMNM/${loci}-mafft-cds.filter/uniques.tre ]; then
                	newline=$(echo ${line} | sed -r 's/busted.fa/uniques_noTree.fa/g' | sed -r 's/busted.tre/uniques.tre/g')
	        else
        	        newline=${line}
	        fi
        	echo ${newline} >> jobs/OUfg_fullReproSet_yesSRV_noMNM.uniques.sh
	fi
done
