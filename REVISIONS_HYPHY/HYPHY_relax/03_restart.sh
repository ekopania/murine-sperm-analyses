#!/bin/bash
#PURPOSE: get loci that did not finish running and restart

#rm done_loci.txt
#ls OUfg_paredSet_noOutgroups/logs/*log | while read file; do 
#	if [ $(grep -c "The following rate distribution was inferred for \*test\* branches" $file) -eq 1 ]; then
#		echo ${file} >> done_loci.txt
#	fi
#done

rm jobs/relax_OUfg_paredSet.REDO.sh
cat jobs/relax_OUfg_paredSet.sh | while read line; do
	if [ $(echo ${line} | cut -c 1) == "#" ]; then
		echo ${line} >> jobs/relax_OUfg_paredSet.REDO.sh
	else
		loci=$(echo ${line} | cut -d "/" -f 13 | cut -d "-" -f 1)
		echo "${loci}"
		if [ $(grep -c "${loci}" done_loci.txt) -eq 0 ]; then
			echo ${line} >> jobs/relax_OUfg_paredSet.REDO.sh
		fi
	fi
done
