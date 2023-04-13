#!/usr/bin/python3

#PURPOSE: If an exon alignment has multiple paralogs for one sample after filtering, keep only the one with the highest exonerate score

import sys, os

sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/GREGG_MURINE_SCRIPTS/");
# Add Gregg's murine lib to the path.

import mseq

#Loop through exonerate output csv: sample-segments-exonerate.csv
#Make a dictionary to get seg-sim values from full sample, exon, and paralog names (should match fasta entries)
#ALL SAMPLES
#exonerate2cds_out = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS-f50/sample-segments-exonerate.csv"
#RTM DATASET
exonerate2cds_out = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/sample-segments-exonerate.csv"
samp_exons = {};
first = True;
for line in open(exonerate2cds_out):
    if first:
        first = False;
        continue;
    line = line.strip().split(',');
    #samp = line[1].strip().split('|')[0];
    #gene = line[1].strip().split("|")[1].split("-")[0];
    #paralog = line[1].strip().split("-")[-1];
    #eid = line[4];
    #pid = line[5];
    # Pretty sure field 9 (seg-sim) is a measure of sequence similarity between the mm10 reference target and the query sequence from the exon capture data; will use this to decide which exon to keep
    seg_sim = line[9];
    #s_e = samp + "_" + eid
    exonerate_entry = line[1];
    samp_exons[exonerate_entry] = { 'seg-sim' : seg_sim};
    #samp_exons[exonerate_entry] = { 'samp-gene' : samp + "|" + gene, 'eid' : eid, 'pid' : pid, 'paralog' : paralog, 'seg-sim' : seg_sim};
    #print(samp, gene, paralog, eid); #For checking/debugging
    #if we haven't run into this sample-exon combination yet, make a new entry in the dictionary
    #if s_e not in samp_exons:
        #samp_exons[s_e] = { 'samp-gene' : samp + "|" + gene, 'eid' : eid, 'pid' : pid, 'paralog' : paralog, 'seg-sim' : seg_sim}; #'paralogs' : []};
    #else: 
        # If this paralog has a better score, keep it and move the old paralog with this sample-exon combination to the "removed" group
        #if seg_sim > samp_exons[s_e]['seg-sim']:
            #print(samp_exons[s_e]['samp-gene'] + "-" + samp_exons[s_e]['paralog']);
            #removed_paralogs.append(samp_exons[s_e]['samp-gene'] + "-" + samp_exons[s_e]['paralog']);
            #samp_exons[s_e] = { 'samp-gene' : samp + "|" + gene, 'eid' : eid, 'pid' : pid, 'paralog' : paralog, 'seg-sim' : seg_sim};
        # Otherwise, put this paralog in the "removed" group and keep the one we already had, because the old one has a better score 
        #else:
            #removed_paralogs.append(line[1]);
        # Either way, we will remove a paralog, so increase the count
        #removed_paralog_count += 1;

#print("Number of paralogs in removed set: " + str(removed_paralog_count));
#print("First few paralogs in removed set:");
#print(removed_paralogs[1:10]);


in_dir = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_byExon_RTM-f16-seq20-site50/nt"
if(os.path.isdir(in_dir + "-removed-paralogs")):
  sys.exit( " * Error: Output directory already exists!");

os.mkdir(in_dir + "-removed-paralogs")

all_fastas = os.listdir(in_dir);

for f in all_fastas:
    #Read in fasta file
    in_fasta = f;
    #print(in_fasta);
    print("Working on fasta: " + in_fasta);
    seqs_orig = mseq.fastaGetDict(in_dir + "/" + in_fasta);
    #print(seqs_orig);
    removed_count=0;
    out_file=in_dir + "-removed-paralogs/" + in_fasta
    #Loop through fasta to figure out which samples have paralogs in this alignment; keep paralog w/ highest exonerate score
    to_write = {};
    paralog_samples=[];
    for title in seqs_orig:
        new_title = title.split(" ")[0];
        new_title = new_title[1:];
        #No mm10 paralogs; just add it to the dict and move on
        if new_title == "mm10":
            short_title = new_title;
            to_write[short_title] = { 't' : title, 'seq' : seqs_orig[title], 'seg-sim' : 100 };
            continue;
        else:
            short_title = new_title[0:new_title.index("|")];
        #Get the exonerat score for this sample, exon, and paralog
        seg_sim = samp_exons[new_title]['seg-sim'];
        #If we haven't seen this sample before, append it to the dictionary
        if short_title not in to_write:
            to_write[short_title] = { 't' : title, 'seq' : seqs_orig[title], 'seg-sim' : seg_sim };
        #If we have seen this sample before but this new paralog has a higher exonerate score, replace dict entry with this paralog
        elif seg_sim > to_write[short_title]['seg-sim']:
            to_write[short_title] = { 't' : title, 'seq' : seqs_orig[title], 'seg-sim' : seg_sim };
            removed_count += 1;
            paralog_samples.append(title);
        #Otherwise, we have seen this sample before and the older version had a better score; keep old version and update removed count
        else:
            removed_count += 1;
            paralog_samples.append(title);
    #For debuggin/checking it worked properly
    print(paralog_samples);
    #Write non-paralgos and paralogs w/ highest score to output fasta in ouput directory
    for s in to_write:
        this_title = to_write[s]['t'];
        mseq.writeSeq(out_file, seqs_orig[this_title], this_title);
    #Print how many paralogs were removed from this fasta
    print("Removed " + str(removed_count) + " samples from fasta");
