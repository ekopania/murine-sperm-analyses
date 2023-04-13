#!/bin/bash
#PURPOSE: Run pare.py to pare tree
#
# Job name:
#SBATCH --job-name=pare_tree
#SBATCH --output=pare_tree-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=32G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro
export PATH=$PATH:/home/ek112884/software/tree-paring/
echo $PATH

#Exempt all branches leading to tips with RTM data
#python /home/ek112884/software/tree-paring/pare.py -t /mnt/beegfs/gt156213e/murinae-seq/docs/data/trees/full_coding_iqtree_astral.cf.rooted.tree -o pared_tree.full_coding_astral.keepRTMtips -i 30 -e rtm_exempt_branches.all_RTM_tips.txt -g 50

#Exempt branches w/ significant optimum shifts based on l1ou
#python /home/ek112884/software/tree-paring/pare.py -t /mnt/beegfs/gt156213e/murinae-seq/docs/data/trees/full_coding_iqtree_astral.cf.rooted.tree -o pared_tree.full_coding_astral.keepOUshifts.tree -i 30 -e rtm_exempt_branches.ou_shifts.txt -g 50

#Exempt all branches leading to tips with RTM data AND branches w/ significant optimum shifts based on l1ou
#python /home/ek112884/software/tree-paring/pare.py -t /mnt/beegfs/gt156213e/murinae-seq/docs/data/trees/full_coding_iqtree_astral.cf.rooted.tree -o pared_tree.full_coding_astral.keepRTMtips_andOUshifts.tree -i 30 -e rtm_exempt_branches.RTMtips_and_OUshifts.txt -g 50

#Start w/ tree that includes only RTM data species; keep branches where there is a shift to lower RTM (i.e., l1ou shifts and their sister clades)
python /home/ek112884/software/tree-paring/pare.py -t /mnt/beegfs/ek112884/murinae/OU_MODELS/pruned_ultrametric_tree.RTM.full_coding_iqtree_astral.cf.rooted.tree -o pared_tree.pruned_ultrametric_RTM.keepOUshifts.tree -i 30 -e rtm_exempt_branches.pruned_ultrametric_RTM.OUshifts.txt -g 0 -b 100 

#Flag explainations; run python /home/ek112884/software/tree-paring/pare.py -h for more info and more options
#  -t TREE_INPUT      Input tree; required
#  -o OUT_DEST        Output directory
#  -g GCF_THRESHOLD   The lower threshold of gene concordance factor to consider for paring at each iteration.
#  -b BL_PERCENTILE   Only branches below this percentile of branch length with be considered for pruning.
#  -i MAX_ITERS       Maximum number of paring iterations 
#  -e EXEMPT          A file containing branches to be exempt from pruning
