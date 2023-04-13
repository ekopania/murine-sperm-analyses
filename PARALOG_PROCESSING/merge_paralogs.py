#!/usr/bin/python3

#NOTES: 9/24/2022
	#If a paralog maps to 2 exons, script is currently appending exon twice - NEED TO FIX THIS; Gregg's script doesn't do this so need to solve that
	#BUT exonerate2cds is currently adding EVERY paralog for a gene to the exonerate2cds output EVEN IF THE PARALOG ONLY MAPS TO EXONS NOT FOUND IN THAT PROTEIN (if found in other proteins of that gene); probably due to modifications I made to Gregg's exonerate2cds script? So need to figure that out
	#UPDATE: If exon appears in multiple proteins, it is only associated with one in the output file; Gregg only included on protein per gene I believe (filtering by txpt w/ most targets); may need to apply this kind of filtering such that only one pid per gene is included BEFORE running exonerate2cds (and possibly earlier in pipeline); some of these proteins are nonsense-mediated decay, etc; this  info is in Ensembl and there are few enough genes that I can filter on this manually

#PURPOSE: Merge "paralogs" that align to different exons (may not be separate paralogs; trying to prevent things from being filtered due to gappiness)

import sys

sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/GREGG_MURINE_SCRIPTS/");
# Add Gregg's murine lib to the path.

import mseq

#Eventually make input directory an argument so that we can do this for different exonerate runs with different parameters, filtering, etc

#Loop through exonerate output csv: sample-segments-exonerate.csv
#exonerate2cds_out = "sample-segments-exonerate.debug_test.csv" #Subset for debugging
#ALL SAMPLES
#exonerate2cds_out = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS-f50/sample-segments-exonerate.csv"
#RTM DATASET
exonerate2cds_out = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/sample-segments-exonerate.csv"
samp_exons = {};
removed_paralog_count = 0;
removed_paralogs = [];
first = True;
for line in open(exonerate2cds_out):
    if first:
        first = False;
        continue;
    line = line.strip().split(',');
    samp = line[1].strip().split('|')[0];
    gene = line[1].strip().split("|")[1].split("-")[0];
    paralog = line[1].strip().split("-")[-1];
    eid = line[4];
    pid = line[5];
    # Pretty sure field 9 (seg-sim) is a measure of sequence similarity between the mm10 reference target and the query sequence from the exon capture data; will use this to decide which exon to keep
    seg_sim = line[9];
    s_e = samp + "_" + eid
    #print(samp, gene, paralog, eid); #For checking/debugging
    #if we haven't run into this sample-exon combination yet, make a new entry in the dictionary
    if s_e not in samp_exons:
        samp_exons[s_e] = { 'samp-gene' : samp + "|" + gene, 'eid' : eid, 'pid' : pid, 'paralog' : paralog, 'seg-sim' : seg_sim}; #'paralogs' : []};
    else: 
        # If this paralog has a better score, keep it and move the old paralog with this sample-exon combination to the "removed" group
        if seg_sim > samp_exons[s_e]['seg-sim']:
            removed_paralogs.append(samp_exons[s_e]['samp-gene'] + "-" + samp_exons[s_e]['paralog']);
            samp_exons[s_e] = { 'samp-gene' : samp + "|" + gene, 'eid' : eid, 'pid' : pid, 'paralog' : paralog, 'seg-sim' : seg_sim};
        # Otherwise, put this paralog in the "removed" group and keep the one we already had, because the old one has a better score 
        else:
            removed_paralogs.append(line[1]);
        # Either way, we will remove a paralog, so increase the count
        removed_paralog_count += 1;

print("Number of paralogs removed: " + str(removed_paralog_count));
print("First few removed paralogs:");
print(removed_paralogs[1:10]);
#Write removed paralogs to file
rp_out = open("removed_paralogs.txt", "w");
for line in removed_paralogs:
    rp_out.write(line);
    rp_out.write("\n");
rp_out.close();

#Some paralogs map to multiple exons - only count these once!

#Loop through dictionary we just made containing sample IDs, exon IDs, gene IDs, and paralog numbers
samp_genes = {};
paralog_tracker = {};
for i in samp_exons:
    s_g = samp_exons[i]['samp-gene'] + "_" + samp_exons[i]['pid'];
    gene = samp_exons[i]['pid'];
    #Use Gregg's fastaGetDict() function from mseq.py to read in the fasta for this gene (protein) as a dictionary
    in_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/nt/" + gene + ".fa";
    #in_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/aa/" + gene + ".fa"
    in_fasta = mseq.fastaGetDict(in_file);
    full_name = samp_exons[i]['samp-gene'] + "-" + samp_exons[i]['paralog'];
    this_key = ">" + full_name;
    #print(in_fasta.keys());
    #Make sure mm10 gets included too
    mm10 = "mm10_" + gene;
    if mm10 not in samp_genes:
        samp_genes[mm10] = { 'pid' : gene, 'seq' : in_fasta['>mm10'] };
    #If this sample-gene combination not already in dictionary, add it to dictionary
    if s_g not in samp_genes:
        samp_genes[s_g] = { 'pid' : gene, 'seq' : "" };
        paralog_tracker[s_g] = { 'paralogs' : [] };
    if samp_exons[i]['paralog'] in paralog_tracker[s_g]['paralogs']:
        print(samp_exons[i]['paralog'] + " already appended; skipping...");
    else:
        #Append sequence to dictionary
        #print("appending" + str(samp_exons[i]));
        samp_genes[s_g]['seq'] = samp_genes[s_g]['seq'] + in_fasta[this_key];
        #print("new samp_genes: " + str(samp_genes[s_g]));
        #Append paralog to paralog tracker
        paralog_tracker[s_g]['paralogs'].append(samp_exons[i]['paralog'])

#Loop through dictionary of sequences and append sequences to the appropriate output fasta file
for i in samp_genes:
    out_file = samp_genes[i]['pid'] + ".removed_paralogs.fa";
    #print(out_file)
    this_title = ">" + i;
    mseq.writeSeq(out_file, samp_genes[i]['seq'], this_title);

print("Done!");
     
#Need to grab only the paralogs that "passed filtering" i.e., didn't have any that matched the same exon and merge ones from same sample
    #OH might be easier to do this by removing problem samples first; then can grab all sequences with the same sample and gene; need to maintain correct exon order though
    #Might need to use Gregg's parseExonerate() function in 08_exonerate_to_cds_2_trimmed.py; can copy the bottom part of this function that writes the file? 
    #Also, eventually want to keep one paralog...best map to mouse? Is there some kind of exonerate score I can go off of?
    #Ok, no, now I'm thinking it'll still be easier to loop through the ones we want to keep; Something like:
        #If there is not already a sequence started for this sample-gene combo, write this sequence from the fasta file and store
        #If there is, append this sequence from the fasta file to the existing seq for this sample-gene combo
        #Probably best to store as dict with gene/protein names as keys; seq for each sample can be a list w/in dict
        #Can then write new fasta for each gene (out of loop, maybe)







##############OLD NOTES#######################


#If this sample, gene, paralog, eid combination exists DO NOT include it; save it to a file of problem genes
    #If doesn't exist, add to a list for this sample/gene combination - MAINTAIN EID ORDER!!!

#Once you've looped through everything for a gene and sample, merge; not sure how I'll know I'm at the end so might need to just wait until eof/out of for loop

#For each unique sample/gene combo
	#Example: Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1
	#If one paralog corresponds to one, different exon
	#Example:
#Ovgp1_ENSMUSE00000381680,Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1-p188,ENSMUSG00000074340,ENSMUST00000000573,ENSMUSE00000381680,ENSMUSP00000000573,252,597,348,78.44827586206897,TRUE,0,0,2
#Ovgp1_ENSMUSE00001220603,Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1-p190,ENSMUSG00000074340,ENSMUST00000000573,ENSMUSE00001220603,ENSMUSP00000000573,413,596,186,91.93548387096774,TRUE,0,0,1
#Ovgp1_ENSMUSE00001246033,Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1-p190,ENSMUSG00000074340,ENSMUST00000000573,ENSMUSE00001246033,ENSMUSP00000000573,687,789,105,100.0,TRUE,0,0,1
#Ovgp1_ENSMUSE00001213690,Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1-p443,ENSMUSG00000074340,ENSMUST00000163626,ENSMUSE00001213690,ENSMUSP00000132424,244,364,123,95.1219512195122,TRUE,0,0,1
#Ovgp1_ENSMUSE00001302293,Abeomelomys_sevia_KUM161018_|_Mus_musculus_Ovgp1-p361,ENSMUSG00000074340,ENSMUST00000163626,ENSMUSE00001302293,ENSMUSP00000132424,311,473,165,96.36363636363636,TRUE,0,0,1 
	#Will probably need to store all sample - gene - paralog # - exon ID info in one set or dict or something
	#If this combo doesn't exist yet, add new one
		#Merge all of these unique paralog/exon combos for the same species/gene into one sequence
		#Need to access the fasta file to do this nt/<prot_id>.fa or maybe aa/<prot_id>.fa
		#IMPORTANT: Gregg seems to have already output exons in the correct order in terms of the protein sequence in his sample-segments-exonerate.csv output from exonerate2cds, so as long as you maintain that order, things should get merged correctly
	#If paralog already has an exon it's associated with for that genes, that's okay; it is spanning multiple exons
		#Include in merge
	#If exon already has a paralog of this gene/sample associated with it, PROBLEM
		#Here is where we likely have paralogs of the same gene/sample overlapping one exon
		#Need to decide which to keep? Or throw out? Is there a way to keep both of them at least until downstream align and filter steps?
		#Lets see how many we actually have in this situation then decide what to do with them
