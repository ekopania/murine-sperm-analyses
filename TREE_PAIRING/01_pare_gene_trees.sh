#!/bin/bash
#PURPOSE: Use newick utils (nw_prune) to prune gene trees to match pared species tree
##Run AFTER paring a tree using pare.py
#
# Job name:
#SBATCH --job-name=prune_gene_trees
#SBATCH --output=prune_gene_tree-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ekopania4@gmail.com
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

dir="pared_tree.pruned_ultrametric_RTM.keepOUshifts.bl100.tree"

echo "Paring gene trees from RTM dataset to match tree in ${dir}"

mkdir "${dir}/GENE_TREES"
#Test with the first few gene trees
#ls /mnt/beegfs/gt156213e/murinae-seq/04-Phylo/reproductive-testes-mass-coding-iqtree/loci/*/*treefile | head | while read file; do
#Run on all RTM dataset gene trees
ls /mnt/beegfs/gt156213e/murinae-seq/04-Phylo/reproductive-testes-mass-coding-iqtree/loci/*/*treefile | while read file; do
	newfile="$(echo ${file} | cut -d "/" -f 9)"
	echo "Working on tree for ${newfile}"
	nw_prune -f "${file}" "${dir}/all-pruned-tips.txt" > "${dir}/GENE_TREES/${newfile}.pared.treefile"
done

echo "Done!"
