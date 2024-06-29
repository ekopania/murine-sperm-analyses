#!/bin/bash
#SBATCH --job-name=getSpermGenes
#SBATCH --output=getSpermGenes-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

echo "#Run HyPhy BUSTED-PH for GO Spermatogenesis genes" > jobs/OUfg_fullReproSet.GOgenesOnly.sh
cat /ix3/nclark/ekopania/MURINAE_REVISIONS/HYPHY_relax/GO_spermatogenesis_proteins.txt | while read loci; do
        line=$(grep ${loci} jobs/OUfg_fullReproSet.sh)
        #Skip empty lines (spermatogenesis protein not in dataset)
        if [ "${line}" == "" ]; then
                continue
        fi
        #If there were duplicate sequences, change the command so we only perform busted-ph on the unique sequences and corresponding tree
        if [ -f OUfg_fullReproSet/${loci}-mafft-cds.filter/uniques.tre ]; then
                newline=$(echo ${line} | sed -r 's/busted-ph.fa/uniques_noTree.fa/g' | sed -r 's/busted-ph.tre/uniques.tre/g')
        else
                newline=${line}
        fi
        echo ${newline} >> jobs/OUfg_fullReproSet.GOgenesOnly.sh
done
