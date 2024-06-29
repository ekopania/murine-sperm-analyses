#!/bin/bash
#SBATCH --job-name=getSpermGenes
#SBATCH --output=getSpermGenes-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

echo "#Run HyPhy RELAX for GO Spermatogenesis genes" > jobs/OUfg_fullReproSet_CPU1.GOgenesOnly.sh
cat GO_spermatogenesis_proteins.txt | while read loci; do
	line=$(grep ${loci} jobs/OUfg_fullReproSet_CPU1.sh)
	#Skip empty lines (spermatogenesis protein not in dataset)
	if [ "${line}" == "" ]; then
		continue
	fi
	#If there were duplicate sequences, change the command so we only perform relax on the unique sequences and corresponding tree
	#Also, fixe the hyphyCPU typo by replacing with hyphy CPU
	if [ -f OUfg_fullReproSet_CPU1/${loci}-mafft-cds.filter/uniques.tre ]; then
		newline=$(echo ${line} | sed -r 's/hyphyCPU/hyphy CPU/' | sed -r 's/relax.fa/uniques_noTree.fa/g' | sed -r 's/relax.tre/uniques.tre/g')
	else
		newline=$(echo ${line} | sed -r 's/hyphyCPU/hyphy CPU/')
	fi
	echo ${newline} >> jobs/OUfg_fullReproSet_CPU1.GOgenesOnly.sh
#	grep ${loci} jobs/OUfg_fullReproSet_CPU1.sh >> jobs/OUfg_fullReproSet_CPU1.GOgenesOnly.sh
done

#CAN'T DO THIS TO WHOLE FILE - some alns have no dups and therefore no uniques.fa files
#Replace "relax.fa" and "relax.tre" with "uniques_noTree.fa" and "uniques.tre")
#sed -i 's/relax.fa/uniques_noTree.fa/g' jobs/OUfg_fullReproSet_CPU1.GOgenesOnly.sh
#sed -i 's/relax.tre/uniques.tre/g' jobs/OUfg_fullReproSet_CPU1.GOgenesOnly.sh
