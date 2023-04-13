#!/usr/bin/python3
############################################################
# For rodent exomes, 04.2020
# Finds the 'best' transcript for each gene based on
# presence in mouse and rat
############################################################

import sys
sys.path.append("/mnt/beegfs/gt156213e/murinae-seq/lib/");
# Add the repo's lib dir to the path.

import os, mcore, gzip, re
from collections import defaultdict

############################################################

mode = "targets";
# Last filter step based on: length or targets

infile = "../repro_exons.protein_coding_only.exon_counts.txt";
gtffile_mouse = "/mnt/beegfs/gt156213e/murinae-seq/Reference-genomes/mm10/Mus_musculus.GRCm38.99.gtf.gz";
gtffile_rat = "/mnt/beegfs/gt156213e/murinae-seq/Reference-genomes/Rnor6/Rattus_norvegicus.Rnor_6.0.99.gtf.gz";

if mode == "targets":
    target_overlaps = "/mnt/beegfs/gt156213e/murinae-seq/Targets/bed/mm10-targets-to-exons-0.9.bed";
    outfilename = "../selected-transcripts-targets.txt"
elif mode == "length":
    outfilename = "../selected-transcripts-length.txt"

#outfilename = "test.txt";

#id_file = "../Reference-genomes/mm10/mm10-ens99-ids.tab";
pid_to_tid = {};
for line in open(infile):
    line = line.strip().split("\t");
    if line[2]:
        pid_to_tid[line[2]] = line[1];
#print(len(pid_to_gid));

overlap_check_file = "/mnt/beegfs/ek112884/murinae/REPRO_PROTEIN_LISTS/prot_list.greenEtal2018.elongating.txt"
overlap_check = open(overlap_check_file, "r").read().split("\n");
subset_len = len(overlap_check);
overlap_check = [ pid_to_tid[pid] for pid in overlap_check if pid in pid_to_tid ];
present_subset_len = len(overlap_check);
master_id_table = "/mnt/beegfs/gt156213e/murinae-seq/02-Annotation-data/mm10-rnor6-master-transcript-id-table.tab"


ds_thresh = 0.5

with open(outfilename, "w") as outfile:
    mcore.runTime("# Rodent exomes -- select mouse trancsripts", outfile);
    mcore.PWS("# Mouse GTF file:        " + gtffile_mouse, outfile);
    mcore.PWS("# Rat GTF file:          " + gtffile_rat, outfile);
    mcore.PWS("# Ensembl ortholog file: " + infile, outfile);
    if mode == "targets":
        mcore.PWS("# Target overlaps file:  " + target_overlaps, outfile);
    mcore.PWS("# Output file:           " + outfilename, outfile);
    mcore.PWS("# --------------", outfile);
    mcore.PWS("# dS threshold:          " + str(ds_thresh), outfile);

    mcore.PWS("# --------------", outfile);

    mcore.PWS("# Checking subset: " + overlap_check_file, outfile);
    mcore.PWS("# Subset genes read : " + str(subset_len), outfile);
    mcore.PWS("# Subset genes found : " + str(present_subset_len), outfile);

    if mode == "targets":
        mcore.PWS("# " + mcore.getDateTime() + " Reading target overlaps...", outfile);
        mouse_transcript_overlaps = {};
        for line in open(target_overlaps):
            line = line.strip().split("\t");
            if line[4] == ".":
                continue;

            gid, tid, eid = line[7].split(";");
            if gid not in mouse_transcript_overlaps:
                mouse_transcript_overlaps[gid] = {};
            if tid not in mouse_transcript_overlaps[gid]:
                mouse_transcript_overlaps[gid][tid] = 0;
            mouse_transcript_overlaps[gid][tid] += 1;

    mcore.PWS("# " + mcore.getDateTime() + " Reading mouse feature lengths...", outfile);
    mouse_transcript_lengths = {};
    exons = {};
    for line in gzip.open(gtffile_mouse):
        line = line.decode();
        if line[0] == "#":
            continue;
        line = line.strip().split("\t");
        feature_type, chrome, start, end, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[8];
        if feature_type == "transcript":
            transcript_len = end - start;
            tid = re.findall('ENSMUST[\d]+', feature_info)[0];
            mouse_transcript_lengths[tid] = { 'chrome' : chrome,  'start' : start, 'end' : end, 'len' : transcript_len, 'num-exons' : 0 };

        if feature_type == "exon":
            mouse_transcript_lengths[tid]['num-exons'] += 1; 

    mcore.PWS("# " + mcore.getDateTime() + " Reading rat feature lengths...");
    rat_transcript_lengths = {};
    rat_pid_to_tid = {};
    for line in gzip.open(gtffile_rat):
        line = line.decode();
        if line[0] == "#":
            continue;
        line = line.strip().split("\t");
        feature_type, chrome, start, end, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[8];
        if feature_type == "transcript":
            transcript_len = end - start;
            tid = re.findall('ENSRNOT[\d]+', feature_info)[0];
            rat_transcript_lengths[tid] = { 'chrome' : chrome, 'start' : start, 'end' : end, 'len' : transcript_len, 'num-exons' : 0 };

        if feature_type == "exon":
            rat_transcript_lengths[tid]['num-exons'] += 1;

        if feature_type == "CDS":
            tid = re.findall('ENSRNOT[\d]+', feature_info)[0];
            pid = re.findall('ENSRNOP[\d]+', feature_info)[0];
            rat_pid_to_tid[pid] = tid;
            
    mcore.PWS("# Mouse transcripts read: " + str(len(mouse_transcript_lengths)), outfile);
    mcore.PWS("# Rat transcripts read:   " + str(len(rat_transcript_lengths)), outfile);
    mcore.PWS("# --------------", outfile);
    
    mcore.PWS("# " + mcore.getDateTime() + " Reading orthologs...", outfile);
    genes = defaultdict(list);
    first = True;
    transcript_count = 0;
    for line in open(infile):
        if first == True:
            orig_headers = line;
            first = False;
            continue;
        if "ENSMUST" in line:
            transcript_count += 1;
        line = line.strip().split("\t");
        gid = line[1];
        genes[gid].append(line[2:]);
    mcore.PWS("# Genes read:      " + str(len(genes)), outfile);
    mcore.PWS("# Transcipts read: " + str(transcript_count), outfile);
    mcore.PWS("# --------------", outfile);

    mcore.PWS("# " + mcore.getDateTime() + " Filtering transcripts...", outfile);
    selected_transcripts = {};

    no_orth = 0;
    no_one2one = 0;
    low_conf = 0;
    ds_filter = 0;
    last_filter = 0;
    no_overlap = 0;
    no_passed_transcripts = 0;

    subset_rm = { 'no-orth' : 0, 'not-one-2-one' : 0, 'low-conf' : 0, 'high-ds' : 0, 'not-targeted' : 0 }

    for gid in genes:
        # print(gid);
        # print(genes[gid]);
        # print("-----")
        if mode == "targets" and gid not in mouse_transcript_overlaps:
            print("No overlap:" + gid);
            no_overlap += 1;
            continue;

        passed_transcripts = [];
        for transcript in genes[gid]:
#            if len(transcript) < 14:
#                no_orth += 1;
#
#                if transcript[0] in overlap_check:
#                    subset_rm['no-orth'] += 1
#                continue;
#
#            if transcript[6] != "ortholog_one2one":
#                no_one2one += 1;
#
#                if transcript[0] in overlap_check:
#                    subset_rm['not-one-2-one'] += 1                
#                continue;
#
#            if transcript[13] != "1":
#                low_conf += 1;
#
#                if transcript[0] in overlap_check:
#                    subset_rm['low-conf'] += 1                
#                continue;
#
#            if transcript[12] == '' or float(transcript[12]) > ds_thresh:
#                ds_filter += 1;
#
#                if transcript[0] in overlap_check:
#                    subset_rm['high-ds'] += 1                
#                continue;

            passed_transcripts.append(transcript);

        if len(passed_transcripts) > 1:

            if mode == "targets":
                max_targets, max_transcript = 0, "";
                for passed_transcript in passed_transcripts:
                    tid = passed_transcript[0];

                    if tid not in mouse_transcript_overlaps[gid]:
                        last_filter += 1;
                        continue;
                    
                    if mouse_transcript_overlaps[gid][tid] > max_targets:
                        max_targets = mouse_transcript_overlaps[gid][tid];
                        max_transcript = passed_transcript;
                    else:
                        if tid in overlap_check:
                            subset_rm['not-targeted'] += 1;
                        last_filter += 1;
                selected_transcripts[gid] = max_transcript;
            # Final selection by target

            if mode == "length":
                max_len, max_transcript = 0, 0;
                for transcript in passed_transcripts:
                # BAD second use of transcript variable within loop... fix later
                    mouse_tid = transcript[0];
                    mouse_len = mouse_transcript_lengths[mouse_tid]['len'];

                    rat_tid = transcript[5];
                    if "ENSRNOP" in rat_tid:
                        rat_tid = rat_pid_to_tid[rat_tid];
                    rat_len = rat_transcript_lengths[rat_tid]['len'];

                    avg_len = (rat_len + mouse_len) / 2;
                    if avg_len > max_len:
                        max_len = avg_len;
                        max_transcript = transcript;
                    else:
                        if mouse_tid in overlap_check:
                            subset_rm['not-targeted'] += 1;
                        last_filter += 1;
                selected_transcripts[gid] = max_transcript;
            # Final selection by length

        elif len(passed_transcripts) == 1:
            max_transcript = passed_transcripts[0];
            selected_transcripts[gid] = max_transcript;

        else:
            no_passed_transcripts += 1;
      

    if mode == "targets":
        mcore.PWS("# Transcripts that do not overlap any targets:                                           " + str(no_overlap), outfile); 
    mcore.PWS("# Transcripts with no rat ortholog:                                                      " + str(no_orth), outfile);
    mcore.PWS("# Transcripts with no one2one rat ortholog:                                              " + str(no_one2one), outfile);
    mcore.PWS("# Transcripts with no high confidence rat ortholog:                                      " + str(low_conf), outfile);
    mcore.PWS("# Transcripts with no rat ortholog below dS threshold:                                   " + str(ds_filter), outfile);
    if mode == "length":
        mcore.PWS("# Transcripts filtered for not being, on average between rat and mouse, longest in gene: " + str(last_filter), outfile);
    if mode == "targets":
        mcore.PWS("# Transcripts filtered for not having max target hits within gene:                       " + str(last_filter), outfile);
    mcore.PWS("# Note: numbers reported are conditional on those transcripts passing the previous filter. Filters reported in order of their coding.");
    mcore.PWS("# --------------", outfile);
    mcore.PWS("# Genes with no passed transcripts:                                                      " + str(no_passed_transcripts), outfile);
    mcore.PWS("# Genes with selected transcript:                                                        " + str(len(selected_transcripts)), outfile);
    mcore.PWS("# --------------", outfile);
    mcore.PWS("# " + mcore.getDateTime() + " Writing selected transcripts to output file...", outfile);

    print(selected_transcripts);
    for t in selected_transcripts:
        outfile.write(selected_transcripts[t][0]);
        outfile.write("\n");
    exit()

    if mode == "length":
        orig_headers = orig_headers.strip() + "\t" + "\t".join(["Mouse chrome", "Mouse start", "Mouse end", "Mouse length", "Mouse exons", "Rat chrome", "Rat start", "Rat end", "Rat length", "Rat exons"]);
    if mode == "targets":
        orig_headers = orig_headers.strip() + "\t" + "\t".join(["Mouse chrome", "Mouse start", "Mouse end", "Mouse length", "Mouse exons", "Target hits", "Rat chrome", "Rat start", "Rat end", "Rat length", "Rat exons"]);
    outfile.write(orig_headers + "\n");
    for gid in selected_transcripts:
        mouse_tid = selected_transcripts[gid][0];
        mouse_gtf = mouse_transcript_lengths[mouse_tid]
        #rat_tid = selected_transcripts[gid][5];
        #if "ENSRNOP" in rat_tid:
        #    rat_tid = rat_pid_to_tid[rat_tid];
        #rat_gtf = rat_transcript_lengths[rat_tid]


        outline = gid + "\t" + "\t".join(selected_transcripts[gid]) + "\t";
        if mode == "length":
            outline += "\t".join([ mouse_gtf['chrome'], str(mouse_gtf['start']), str(mouse_gtf['end']), str(mouse_gtf['len']), str(mouse_gtf['num-exons']) ]) + "\t";
        if mode == "targets":
            outline += "\t".join([ mouse_gtf['chrome'], str(mouse_gtf['start']), str(mouse_gtf['end']), str(mouse_gtf['len']), str(mouse_gtf['num-exons']), str(mouse_transcript_overlaps[gid][mouse_tid]) ]) + "\t";
        #outline += "\t".join([ rat_gtf['chrome'], str(rat_gtf['start']), str(rat_gtf['end']), str(rat_gtf['len']), str(mouse_gtf['num-exons']) ]);
        outline += "\t".join([ str(mouse_gtf['num-exons']) ]);
        outfile.write(outline + "\n");

print(subset_rm);
