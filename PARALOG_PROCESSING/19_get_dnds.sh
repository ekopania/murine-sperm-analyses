#!/bin/bash
#PURPOSE: Get dN/dS value from M0 run
#SBATCH --job-name=get_m0_dnds
#SBATCH --output=get_m0_dnds-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2M
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
# Commands to run:

myPath="/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m0"

echo "protID	dN.dS" > "${myPath}/dnds.RTM_byExon-f16-seq20-site50_m0.txt"
ls -d ${myPath}/ENS* | while read dir; do
	#echo "${dir}"
	name=$(echo "${dir}" | cut -d "/" -f 11 | cut -d "-" -f 1)
	dnds=$(grep "omega (dN/dS)" "${dir}/codeml.out" | awk '{print $4}')
	if [ -z ${dnds} ]; then
		dnds="NA"
	fi
	echo "${name}	${dnds}" >> "${myPath}/dnds.RTM_byExon-f16-seq20-site50_m0.txt" 
done

echo "Done!"
