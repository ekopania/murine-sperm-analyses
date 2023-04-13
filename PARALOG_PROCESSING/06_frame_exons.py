#!/usr/bin/python3
###########################################################
# For rodent exomes, 02.2021
# Takes CDS exons and trims off split codons at splice sites
###########################################################

import sys
sys.path.append("/home/ek112884/software/core/corelib/");
# Add Gregg's corelib dir to the path.

import os, core, math, argparse, re, subprocess, multiprocessing as mp
from Bio.Seq import Seq

###########################################################
# Globals

ref = "/mnt/beegfs/ek112884/REFERENCE_DIR/mm10.fa";

###########################################################
# Functions

def sortCoords(coord_list, strand, exclude=[]):
# Takes a set of coordinates and sorts them in ascending order.
# Returns a list of the indices of the sorted entries in the original
# list, which is useful when multiple lists need to be sorted in the
# same way (i.e. start and end coordinates of genomic fragments).
    exclude = list(set(exclude));
    reverse = False;
    if strand == "-":
        reverse = True;


    sorted_coords = sorted(coord_list, reverse=reverse);
    sorted_inds = [];

    for c in sorted_coords:
        start_ind = coord_list.index(c);
        if start_ind in exclude:
            exclude.remove(start_ind);
        if start_ind not in exclude:
            sorted_inds.append(coord_list.index(c));

    return(sorted_inds);

################################

def runFaidx(c, s, e, strand):
# Given a set of genomic coordinates, runs faidx to retrieve that sequence.
    cmd = "samtools faidx " + ref + " " + c + ":" + s + "-" + e;
    if strand == "-":
        cmd += " -i";     
    #print(cmd);
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    assert cmd_result.stderr.decode() == "", "\nsamtools returned an error with the following command:\n" + cmd + "\n\n " + cmd_result.stderr.decode() + "\n";

    samtools_raw_out = cmd_result.stdout.decode().split("\n");
    header = samtools_raw_out[0];
    seq = "".join(samtools_raw_out[1:]);

    return header, seq;

################################

def writeSeq(header, seq, outfilename):
# Writes a sequence to a file in FASTA format.
    with open(outfilename, "w") as outfile:
        outfile.write(header + "\n");
        outfile.write(seq + "\n");

################################

def getExons(pid_list, pid_dict, out_dir_nt, out_dir_aa, debug=False):
    ##########
    # Begin protein block
    
    num_prots = 0;
    stop_codon_prots = [];
    for pid in pid_list:
    # Retrieve alignments for every protein in the dataset.

        if debug:
            print(pid);
        # Debug statements.   

        cur_chr, cur_strand = pid_dict[pid]['chrome'], pid_dict[pid]['strand'];
        gid, tid = pid_dict[pid]['gid'], pid_dict[pid]['tid'];
        # Get info for current protein.

        pid_header = ">" + "|".join([gid, tid, pid]);
        # Make the header for the protein.

        sorted_indices = sortCoords(pid_dict[pid]['estarts'], cur_strand);
        # Sort the exons depending on strand.

        first_exon = True;
        last_exon = False;
        running_len = 0;
        prev_trim = 0;
        next_trim = 0;
        # These variables keep track of things between exons.
        for i in sorted_indices:
        # Go through every exon in the protein.
            
            if i == sorted_indices[-1]:
                last_exon = True;

            cur_id, cur_start, cur_end = pid_dict[pid]['eids'][i], str(pid_dict[pid]['estarts'][i]), str(pid_dict[pid]['eends'][i]);
            # Get info for the current exon.

            # print(cur_id, cur_start, cur_end);
            # cur_len = int(cur_end) - int(cur_start);
            # cur_target = "-".join([pid_dict[pid]['chrome'], cur_start, cur_end]);
            
            # r = re.compile(cur_id);
            # file_matches = list(filter(r.match, seq_files));
            # print(file_matches);
            # assert len(file_matches) == 1, "\nNO MATCHING FILES FOUND: " + cur_id;


            # cur_seq_file = os.path.join(indir, file_matches[0]);
            # print(cur_seq_file);
            # assert os.path.isfile(cur_seq_file), "\nCANNOT FIND SEQ FILE: " + cur_seq_file;
            
            # cur_header, cur_seq = core.fastaGetLists(cur_seq_file);
            # cur_header = cur_header[0];
            # cur_seq = cur_seq[0];

            # assert len(cur_seq) == cur_len, "\nUNEXPECTED SEQ LEN: " + cur_id + "\n\tEXPECTED: " + str(cur_len) + "\n\tACTUAL: " + str(len(cur_seq));
            cur_header, cur_seq = runFaidx(cur_chr, cur_start, cur_end, cur_strand);
            # Get the exon sequence.

            if debug:
                print("\t" + cur_id);
                print("\t\tHEADER         :", cur_header);
                print("\t\tORIG SEQ       :", cur_seq);
                print("\t\tORIG SEQLEN    :", len(cur_seq));
                print("\t\tPREV BACK TRIM :", prev_trim);
                print("\t\tTO TRIM FRONT  :", next_trim);
            # Debug statements.
            
            cur_seq = cur_seq[next_trim:];
            # Trim the sequence based on the amount trimmed in the last exon.

            if debug:
                print("\t\tINT SEQ:        ", cur_seq);
            # Debug statements.

            seqlen = len(cur_seq);
            modlen = seqlen % 3;
            # Get the intermediate trimmed sequence length.

            if debug:
                print("\t\tINT SEQLEN     :", seqlen);
                print("\t\tTO TRIM BACK   :", modlen);
            # Debug statements.

            if modlen != 0:
                cur_seq = cur_seq[:-modlen];
                next_trim = 3 - modlen;
            else:
                next_trim = 0;
            prev_trim = modlen;
            # Check if the intermediate length needs to be trimmed at the end.
                   
            if debug:
                print("\t\tNEXT FRONT TRIM:", next_trim);
            # Debug statements.

            final_len = len(cur_seq);
            if final_len % 3 != 0:
                print("OUT OF FRAME EXON!", pid)
                return "exit";
            # Get the final exon length. If out of frame, report and exit.

            running_len += final_len;
            # Add exon length to full protein length.

            if debug:
                print("\t\tFINAL SEQ      :", cur_seq);
                print("\t\tFINAL SEQLEN   :", final_len);
            # Debug statements.

            cur_seq_aa = str(Seq(cur_seq).translate());
            # Translate the sequence for AA output.

            if debug:
                print("\t\tAA SEQ         :", cur_seq_aa);
            # Debug statements.

            if len(cur_seq) > 3 and last_exon and cur_seq_aa[-1] == "*":
                if debug:
                    print("\t\tFINAL POS STOP CODON... REMOVING")
                cur_seq = cur_seq[:-3];
                cur_seq_aa = cur_seq_aa[:-1];
            # CHECK LAST POS FOR STOP HERE #####

            final_header = pid_header + "|" + cur_id + "|" + cur_header.replace(">", "");
            # Add exon info to output header

            if debug:
                print("\t\tOUT HEADER     :", final_header);
            # Debug statements.

            nt_out = os.path.join(out_dir_nt, cur_id + "-nt-trimmed.fa");
            writeSeq(final_header, cur_seq, nt_out);
            aa_out = os.path.join(out_dir_aa, cur_id + "-aa-trimmed.fa");
            writeSeq(final_header, cur_seq_aa, aa_out);
            # Write both NT and AA sequences.

            if "*" in cur_seq_aa:
                print("STOP CODON!", pid);
                stop_codon_prots.append([pid, cur_id]);
                #return "exit";
            # Exit if there is a stop codon.

        if running_len % 3 != 0:
            print("OUT OF FRAME CDS!", pid)
            return "exit";
        # Exit if the protein is out of frame.

        num_prots += 1;
        # if num_prots > 1000:
        #     return stop_codon_prots;
        # Exit early for the test set.

    return stop_codon_prots;

###########################################################
# Main block

#annotation_file = "../02-Annotation-data/Mus-selected-sequences_metadata_samp-counts_ratids-TESTSET.csv";
#annotation_file = "../02-Annotation-data/Mus-selected-sequences_metadata_samp-counts_ratids.csv";

annotation_file = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/repro_exons_full.txt"
#annotation_file = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/repro_exons_TEST.txt"
# Table matching protein, transcript, gene, and exon IDs

#indir = "../Targets/seq/mm10-target-exon-overlaps/";
outdir_nt = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/repro_genes-selected-cds-nt-trimmed/";
outdir_aa = "/mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/repro_genes-selected-cds-aa-trimmed/";
# Input and output directories.

#seqfiles = os.listdir(indir);

logfilename = os.path.join("logs", "frame-exons.log");
# Log file for the run.

pad = 25
procs = 1;
# Job vars.

with open(logfilename, "w") as logfile:
    core.PWS("# Reading annotation file: " + annotation_file, logfile);
    mus_pids = {};
    # The mus_pids dict contains the conversion between mouse IDs and target/contig IDs.

    num_targets = 0;
    first = True;
    for line in open(annotation_file):
        if first:
            first = False;
            continue;
        line = line.replace("\"", "").strip().split("\t");
        targ, egn, gid, tid, pid, eid = line[0].replace("\"", "")+"_"+line[4].replace("\"", ""), line[0].replace("\"", ""), line[1].replace("\"", ""), line[2].replace("\"", ""), line[3].replace("\"", ""), line[4].replace("\"", "");
        #targ, gid, tid, pid, eid = line[0].replace("\"", ""), line[5].replace("\"", ""), line[9].replace("\"", ""), line[13].replace("\"", ""), line[17].replace("\"", "");
        if pid == "NA":
            continue;
        exon_start, exon_end = int(line[10]), int(line[11]);
        chrome, strand = "chr" + line[5], line[12];
        if pid not in mus_pids:
            mus_pids[pid] = { 'gid' : gid, 'tid' : tid, 'eids' : [], 'estarts' : [], 'eends' : [], 'targets' : [], 'chrome' : chrome, 'strand' : strand };

        mus_pids[pid]['eids'].append(eid);
        mus_pids[pid]['estarts'].append(exon_start);
        mus_pids[pid]['eends'].append(exon_end);
        mus_pids[pid]['targets'].append(targ);
        num_targets += 1;
    core.PWS("# Protein IDs read: " + str(len(mus_pids)), logfile);
    core.PWS("# Total targets: " + str(num_targets));
    core.PWS("# ----------------", logfile);
    # Reads the mouse and target annotation info.

    ##########################
    # Main parsing block -- splits out into multiple processes.

    # written_main, skipped_main, no_hits_main = 0,0,[];
    # tracker_main = {};
    # The main count and tracker variables. The tracker will be updated for each chunk of proteins.

    pids = list(mus_pids.keys());
    pids_per_proc = int(math.ceil(len(pids) / procs));
    pid_lists = list(core.chunks(pids, pids_per_proc));
    # Splitting the protein IDs into even chunks for parallelization

    pool = mp.Pool(processes=procs);
    # Starting the pool of processes

    stop_codon_prots_main = [];    

    for result in pool.starmap(getExons, ( (pl, mus_pids, outdir_nt, outdir_aa, True) for pl in pid_lists )):
        if result == "exit":
            sys.exit();
        stop_codon_prots_main += result;
        # tracker_main.update(result[0]);
        # written_main += result[1];
        # skipped_main += result[2];
        # no_hits_main += result[3];
    # Pass each chunk of IDs to the parse function and collect variables
    core.PWS(core.spacedOut("\n# Exons with stop codons:", pad) + str(len(stop_codon_prots_main)), logfile);
    with open("mm10-exons-w-stops.txt", "w") as outfile:
        for ids in stop_codon_prots_main:
            outfile.write(",".join(ids) + "\n");
    #core.PWS(core.spacedOut("# Files written:", pad) + str(written_main), logfile);
    #core.PWS(core.spacedOut("# Files skipped:", pad) + str(skipped_main), logfile);
    core.PWS("# ----------------", logfile);
    ##########################
