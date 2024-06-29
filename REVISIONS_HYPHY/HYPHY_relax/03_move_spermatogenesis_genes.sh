#!/bin/bash
#PURPOSE: Move spermatogenesis loci to a separate directory

mkdir GOspermatogenesisProteins_OUfg_fullReproSet/
ls OUfg_fullReproSet_CPU1/*json | while read file; do
	echo ${file}
	dir=$(echo ${file} | cut -d "." -f 1-2)
	echo ${dir}
	mv ${file} GOspermatogenesisProteins_OUfg_fullReproSet/
	mv ${dir} GOspermatogenesisProteins_OUfg_fullReproSet/
done

mkdir GOspermatogenesisProteins_OUfg_fullReproSet/logs
mv OUfg_fullReproSet_CPU1/logs/*log GOspermatogenesisProteins_OUfg_fullReproSet/logs
