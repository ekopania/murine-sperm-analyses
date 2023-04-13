#!/bin/bash
##PURPOSE: use a LRT to compare lnL values of bs test and bs null, and of bs null vs M1a; use the chi2 function of PAML to test for significant LRT result (pos selection; df=2)
#SBATCH --job-name=bsTest_LRT
#SBATCH --output=bsTest_LRT-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=ekopania4@gmail.com
#SBATCH --partition=good_lab_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Partition:
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
# Commands to run:

module load paml
module load R/4.1.1

echo "protID	bs_lnL	bsNull_lnL	m1a_lnL	bs_LRT	bs_p	relaxed_LRT	relaxed_p" > /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.bsTest.top10.txt
ls -d /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/RTM_SET_pared_branchSite_reproGenes_top10/* | while read dir;
do
	gene=$(echo "${dir}" | cut -d "/" -f 11 | cut -d "-" -f 1)
	echo "${gene}"
	lnL=$(grep "lnL" "${dir}/codeml.out" | cut -d ":" -f 4 | awk '{print $1}')
	echo "${lnL}"
	lnL_null=$(grep "lnL" /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/RTM_SET_pared_branchSiteNull_reproGenes_top10/${gene}-mafft-cds.filter/codeml.out | cut -d ":" -f 4 | awk '{print $1}')
	echo "${lnL_null}"
	#Get the difference and also replace e with ^ if the answer is output in sci not so that bc can read it
	LRT=$(awk 'function abs(v) {return v < 0 ? -v : v} {print 2*($2-$1)}' <<< "$lnL_null $lnL")
	r=$(awk -v a="-2e-06" -v b="0" 'BEGIN{print (a<b)?1:0}')
	if [ $(awk -v a="${LRT}" -v b="0" 'BEGIN{print (a<b)?1:0}') -eq 1 ]; then
		LRT="NA"
	fi
	if [ -z "${lnL_null}" ]; then #check if length zero
		LRT="NA"
	fi
	p=$(chi2 2 "${LRT}" | awk '{print $8}')
	echo "${p}"
	lnL_m1a=$(grep "lnL" /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/RTM_SET_pared_M1a/${gene}-mafft-cds.filter/codeml.out | cut -d ":" -f 4 | awk '{print $1}')
        echo "${lnL_null}"
	LRT_relaxed=$(awk 'function abs(v) {return v < 0 ? -v : v} {print 2*($2-$1)}' <<< "$lnL_m1a $lnL_null")
	if [ $(awk -v a="${LRT_relaxed}" -v b="0" 'BEGIN{print (a<b)?1:0}') -eq 1 ]; then
                LRT_relaxed="NA"
        fi
        if [ -z "${lnL_null}" ]; then
                LRT_relaxed="NA"
        fi
        echo "${LRT_relaxed}"
	p_relaxed=$(chi2 2 "${LRT_relaxed}" | awk '{print $8}')
	echo "${p_relaxed}"
	printf "%s\t%f\t%f\t%f\t%s\t%f\t%s\t%f\n" "${gene}" "${lnL}" "${lnL_null}" "${lnL_m1a}" "${LRT}" "${p}" "${LRT_relaxed}" "${p_relaxed}" >> /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/DNDS_BY_SPERM_STAGE/PAML/pos_selec.bsTest.top10.txt
done

Rscript 16_multiple_test_correct.bsTest.r

echo "Done!"
