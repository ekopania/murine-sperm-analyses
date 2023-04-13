#!/usr/bin/python
############################################################
# 09.2022
# Filter paralog FASTA files based on BLAST results.
# BLAST results must be in -outfmt 6 
###########################################################

import sys
sys.path.append("/home/ek112884/software/core/corelib/");
# Add Gregg's corelib dir to the path.

import os, argparse, core
from collections import defaultdict
from Bio import SeqIO

############################################################

# IO options
parser = argparse.ArgumentParser(description="Exonerate command generator");
parser.add_argument("-i", dest="in_fasta", help="The input fasta with all paralog sequences corresponding to a given gene.", default=False);
#parser.add_argument("-s", dest="score_mode", help="The score type to maximize by. One of: len, eval, bit", default=False);
parser.add_argument("-e", dest="eval_cutoff", help="The max evalue to consider a hit. Default: 1e-3", type=float, default=1e-3);
parser.add_argument("-o", dest="out_fasta", help="Desired output fasta file name.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
args = parser.parse_args();

#Input and error checkiing
if args.in_fasta == False:
    parser.print_help();
    sys.exit(" * ERROR 1: Must speficy an input fasta file with -i. See help menu with -h.");
if args.out_fasta == False:
    parser.print_help();
    sys.exit(" * ERROR 2: Must speficy an output fasta file with -o. See help menu with -h.");
if os.path.isfile(args.out_fasta) and not args.overwrite:
    sys.exit(" * ERROR 3: Output file already exists! Explicity specify --overwrite to overwrite it.");

#Get external gene name based on input file
gene_name = (args.in_fasta.strip().split("/")[-1]).split(".")[0]
print("Working on gene: " + gene_name)

#Read in exon to gene information
exon_info = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/repro_exons_full.txt"

#Get all exons associated with gene we are working on
gene_exons = []
first = True;
for line in open(exon_info):
    if first:
        first = False;
        continue;
    line = line.strip().split("\t");
    this_gene = line[0];
    if this_gene == gene_name:
        if(line[4] not in gene_exons):
            gene_exons.append(line[4]);
    
print("Number of exons associated with this gene: " + str(len(gene_exons)))
print(gene_exons)

#Loop through every BLAST output for every exon associated with this gene and keep paralog sequences that pass BLAST filtering for any of these exons
print("Keeping paralogs with BLAST e-value < " + str(args.eval_cutoff))
keep = []
for e in gene_exons:
    blast_tab = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/BLAST_RESULTS/" + e + "-nt-trimmed.blastn.to" + gene_name + ".out"
    for line in open(blast_tab):
        line = line.strip().split("\t");
        if float(line[10]) < args.eval_cutoff:
            if line[1] not in keep:
                keep.append(line[1]);

print("Number of paralog sequences in the gene fasta that passed BLAST filtering:" + str(len(keep)))
#print(keep)

#Filter fasta file
keep_set = set(keep); #Biopython manual says this is faster as a set
records = (r for r in SeqIO.parse(args.in_fasta, "fasta") if r.id in keep_set);
count = SeqIO.write(records, args.out_fasta, "fasta");
print("Saved %i records from %s to %s" % (count, args.in_fasta, args.out_fasta));
if count < len(keep):
    print("Warning %i IDs not found in %s" % (len(keep) - count, args.in_fasta));
    

