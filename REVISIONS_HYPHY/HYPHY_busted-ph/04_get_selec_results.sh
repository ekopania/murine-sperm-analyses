#!/bin/bash
#PURPOSE: Get loci with evidence for positive selection associated with phenotype/trait
#SBATCH --job-name=get_selec
#SBATCH --output=get_selec-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M


#GO spermatogenesis proteins; selection associated with trait ONLY
#echo "protID	selec_result" > GOspermatogenesisProteins_OUfg_fullReproSet.selecAssocTrait.txt
#grep -w "Selection is associated with the phenotype / trait$" GOspermatogenesisProteins_OUfg_fullReproSet/logs/*log | while read line; do
#	prot=$(echo ${line} | cut -d "/" -f 3 | cut -d "-" -f 1)
#	res=$(echo ${line} | cut -d ":" -f 2)
#	echo "${prot}	${res}" >> GOspermatogenesisProteins_OUfg_fullReproSet.selecAssocTrait.txt
#done

#ALL proteins; ALL results
#echo "protID    selec_result" > OUfg_fullReproSet.selecAssocTrait.txt
#ls GOspermatogenesisProteins_OUfg_fullReproSet/logs/*log | while read file; do 
#	prot=$(echo "${file}" | cut -d "/" -f 3 | cut -d "-" -f 1)
#	res=$(tail -1 ${file})
#	echo "${prot}	${res}" >> OUfg_fullReproSet.selecAssocTrait.txt
#done
#ls OUfg_fullReproSet/logs/*log | while read file; do
#        prot=$(echo "${file}" | cut -d "/" -f 3 | cut -d "-" -f 1)
#        res=$(tail -1 ${file})
#        echo "${prot}   ${res}" >> OUfg_fullReproSet.selecAssocTrait.txt
#done

#Get error loci (should just be the spermatogenesis genes that I moved)
grep "error" OUfg_fullReproSet.selecAssocTrait.txt | awk '{print $1}' > error_loci.txt
#echo "#Rerun loci that gave errors the first time" > jobs/OUfg_fullReproSet.REDO.sh
#cat error_loci.txt | while read line; do
#	grep "${line}" jobs/OUfg_fullReproSet.sh >> jobs/OUfg_fullReproSet.REDO.sh
#done

#Get loci WITHOUT errors
grep -v "error" OUfg_fullReproSet.selecAssocTrait.txt > OUfg_fullReproSet.selecAssocTrait.noErrors.txt
