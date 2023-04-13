#!/bin/bash
#SBATCH --job-name=ssm_perm_logRTM_binary
#SBATCH --output=ssm_perm_logRTM_binary-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.kopania@genetics.utah.edu
#SBATCH --account=clarkn-rw
#SBATCH --partition=clarkn-shared-rw
##SBATCH --account=redwood-shared-short
##SBATCH --partition=redwood-shared-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=7-0

#source ~/software/anaconda/anaconda3/bin/activate
#conda activate r4

#module load libxml2
module load R/4.1.1

echo "Modules loaded; starting command"

#Rscript /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/09_permulations_only.binary.SSM.r
Rscript /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/09_permulations_only.binary.SSM.r 20 
#Rscript /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/12.3_permulationsWenrichment.binary.forParallel.R 1
#parallel -j 40 < /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/RERconverge/NEW_RTM_DATA_NOV2022/jobs/ssm_perm_logRTM_binary.sh
