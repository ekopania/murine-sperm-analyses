#!/bin/bash
##PURPOSE: use a LRT to compare lnL values of M1a and M2a, use the chi2 function of PAML to test for significant LRT result (pos selection; df=2)
#SBATCH --job-name=m1aVm2a_LRT
#SBATCH --output=m1aVm2a_LRT-%j.out
##SBATCH --mail-type=ALL
###SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
# Commands to run:

module load paml

rm /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.m1aVm2a.txt
ls -d /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m2a/* | while read dir;
do
	gene=$(echo "${dir}" | cut -d "/" -f 11 | cut -d "-" -f 1)
	echo "${gene}"
	lnL=$(grep "lnL" "${dir}/codeml.out" | cut -d ":" -f 4 | awk '{print $1}')
	echo "${lnL}"
	lnL_null=$(grep "lnL" /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/RTM_byExon-f16-seq20-site50_m1a/${gene}-mafft-cds.filter/codeml.out | cut -d ":" -f 4 | awk '{print $1}')
	echo "${lnL_null}"
	LRT=$(awk 'function abs(v) {return v < 0 ? -v : v} {print 2*($2-$1)}' <<< "$lnL_null $lnL")
	if [ $(echo "${LRT} < 0" | bc) -eq 1 ]; then
		LRT="NA"
	fi
	if [ -z "${lnL_null}" ]; then
		LRT="NA"
	fi
	echo "${LRT}"
	p=$(chi2 2 "${LRT}" | awk '{print $8}')
	echo "${p}"
	printf "%s\t%f\t%f\t%f\t%f\n" "${gene}" "${lnL}" "${lnL_null}" "${LRT}" "${p}" >> /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/SELEC_TESTS/pos_selec.m1aVm2a.txt
done

Rscript 13_multiple_test_correct.m1aVm2a.r

echo "Done!"
