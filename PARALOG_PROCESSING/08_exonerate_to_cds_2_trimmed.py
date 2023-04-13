#!/usr/bin/python3
###########################################################
# For rodent exomes, 12.2020
# Modified by E. Kopania 08.2022
# Extracts and concatenates aligned segments from 
# exonerate .gff output
#
# Definitions:
# target - an exome contig sequence predicted for a given protein.
# query  - the mouse protein sequence.
###########################################################

import sys

sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/software/core/corelib/");
# Add Gregg's corelib dir to the path.
sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/GREGG_MURINE_SCRIPTS/");
# Add Gregg's murine lib to the path.

import os, core, math, argparse, subprocess, multiprocessing as mp, mseq
from Bio.Seq import Seq
from difflib import SequenceMatcher
from glob import glob

###########################################################
# Globals
aln_type_str = ["match", "start to match", "match to end", "delete", "match to delete", "delete to match", "insert", "match to insert", "insert to match",
                "frameshift open 1 p2d", "frameshift open 2 p2d", "frameshift close 3 p2d", "frameshift close 0 p2d",
                "frameshift open 2 target", "frameshift close 0 target", "frameshift open 1 target",
                "frameshift open 1 query", "frameshift close 0 query",
                "frameshift close 3 target", 
                "frameshift open 2 query", "frameshift close 3 query",
                "(START) to intron 0:0 phase-T", "target intron loop 0:0 phase-T", "intron 0:0 phase-T to (END)",
                "(START) to phase1pre phase-T", "(START) to intron 1:2 phase-T", "target intron loop 1:2 phase-T", "intron 1:2 phase-T to (END)", "phase1post phase-T to (END)",                
                "(START) to phase2pre phase-T", "(START) to intron 2:1 phase-T", "target intron loop 2:1 phase-T", "intron 2:1 phase-T to (END)", "phase2post phase-T to (END)"]
# List of valid exonerate alignment names. Add as encountered.

intron_start_str = ["(START) to intron 0:0 phase-T", "(START) to intron 1:2 phase-T", "(START) to intron 2:1 phase-T"];
intron_split_str = ["(START) to phase1pre phase-T", "(START) to phase2pre phase-T"];
intron_end_str = ["phase1post phase-T to (END)", "phase2post phase-T to (END)"];
# List of valid exonerate intron alignment starts and ends for split codons.

###########################################################
# Functions

def readSeqInfo(seqfile, tabfile):
# This function takes a multi-sequence FASTA file and a tab file with annotation info and reads
# the sequences and whether or not they are expected. The tab file should be from gxf_parse.py 
# which is my GFF/GTF parsing script.

    expecteds = {};
    # Exons that are found in the GFF file and their coding exons.

    first = True;
    for line in open(tabfile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        # Skip the first line and comment lines.

        line = line.strip().split("\t");
        ftype, tid, start, end, num_cds = line[1], line[6].split(";")[1], int(line[3])-1, int(line[4]), line[7];
        # Split the line and unpack info.

        if ftype != "protein_coding exon":
            continue;
        # If the exon is not protein coding, skip.

        if tid not in expecteds:
            expecteds[tid] = { 'num-cds' : 0, 'tot-len' : 0, 'starts' : [], 'ends' : [] };
        # If the transcript hasn't been added to the expected dict, initialize here.

        expecteds[tid]['num-cds'] += 1;
        cds_len = end - start;
        expecteds[tid]['tot-len'] += cds_len;
        expecteds[tid]['starts'].append(start);
        expecteds[tid]['ends'].append(end);
        # Save info for the current exon.

    core.PWS("# Protein coding transcripts read: " + str(len(expecteds)), logfile);
    # Read tab file.

    # for tid in expecteds:
    #     start = min(expecteds[tid]['starts']);
    #     end = max(expecteds[tid]['ends']);
    #     tlen = end - start;
    #     expecteds[tid]['len'] = tlen;
    # # Get CDS lengths.

    seqs_orig = mseq.fastaGetDict(seqfile);
    # Read the sequences into a dictionary.
    seqs = {};
    for title in seqs_orig:
        new_title = title.split(" ")[0];
        new_title = new_title[1:new_title.index(".")];
        if new_title not in expecteds:
            continue;
        seqs[new_title] = seqs_orig[title];
    # Check if each sequence has an expected transcript ID. If not, skip. If so, add it to final dict of sequences.
    core.PWS("# Coding sequences read: " + str(len(seqs)), logfile);
    # Read sequences;

    for tid in seqs:
        assert abs(expecteds[tid]['tot-len'] - len(seqs[tid])) <= 10, "\nUnexpected sequence length:\n" + tid + "\nExpected: " + str(expecteds[tid]['tot-len']) + "\nObserved: " + str(len(seqs[tid])) + "\n" + str(expecteds[tid]);
    # Check that the lengths match.

    return seqs, expecteds;
    # Return.

################################

def splitSeq(seq, ind_list, seq_type):
# This function takes a sequence (amino acid or codon) and a list of indices
# and splits that sequence into a list at those indices. Residues/codons at 
# those indices are removed.

    split_seqs = [];
    # The list of split sequences to be determined and returned.

    if seq_type in ["aa", "nt"]:
        interval = 1;
    elif seq_type == "codon":
        interval = 3;
    # The interval to increment the index depends on the sequence type.

    seq_chunks, prev_ind, first_chunk = [], 0, True;
    for ind in ind_list:
        seq_chunks.append(seq[prev_ind:ind]);
        prev_ind = ind + interval;
    # For each index, splice out the string and append to a list.

    if prev_ind < len(seq):
        seq_chunks.append(seq[prev_ind:]);
    # Add on the last segment of the sequence.

    return seq_chunks;

################################

def percIdent(seq1, seq2):
# This function calculates the percent of positions in two aligned sequences
# that are identical.

    # print(seq1);
    # print(seq2);

    ident_pos = 0;
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            ident_pos += 1;
    
    return ident_pos / len(seq1);

################################

def codonComplete(codon):
    if codon == "-":
        codon = "---";
    # For gaps in the target sequence, expand to triple gaps for the codon.

    while len(codon) != 3:
        codon += "N";
    # For frameshifts, insert Ns to complete the codon.


    return codon;

################################

def parseExonerate(pid_list, mus_pids, indir, outdir, filter_thresh, noseq, debug=True):
# This function parses the exonerate output files to retrieve sequences and alignments.
# 1. Gets all predicted segments for each target.
#    - Alignments read to incorporate frameshifts (insert Ns) and introns (concatenate split codons).
# 2. Chooses longest predicted segment.
# 3. Checks chosen segement for stop codons.
# 4. If stop codons are found, splits sequence by stop codon and chooses chunk with highest
#    similarity to mouse sequence.
# 5. Concatenates all chosen segements from a given protein and writes to file.

    target_tracker = {};
    # Info for every target.

    sample_segment_tracker = {};
    # Info for every sample.

    sample_target_tracker = {};
    # Info for every target in every sample.

    written, skipped = 0, 0;
    # Output counters.

    debug_condition, debug_outcome = True, False;
    print("Debug", debug);

    ##########
    # Begin protein block
    for pid in pid_list:
    # Retrieve alignments for every protein in the dataset.

        mtid = mus_pids[pid]['tid'];
        # if mtid not in mexp:
        #     skipped += 1;
        #     continue;
        #rpid = mus_pid_to_rat_pid[pid];
        #rtid = rat_pid_to_tid[rpid];
        # Get the corresponding rat and mouse reference IDs.

        sample_target_tracker[pid] = {};

        sequence_tracker = {};
        sequence_tracker["mm10"] = {'nt-seq' : "", 'aa-seq' : ""}
        # Tracker dicts for the current protein.
        # target_tracker is the main dict that tracks targets.
        # sample_tracker is the protein specific tracker that tracks concatenated targets for the current protein.

        # if debug:
        #     print(pid);
        #     print(mus_pids[pid]['targets'])

        ##########
        # Begin target block
        for target in mus_pids[pid]['targets']:
        # Go through every target in the current protein.

            target_index = mus_pids[pid]['targets'].index(target);
            # Get the index of the current target in the current protein

            meid = mus_pids[pid]['eids'][target_index];
            # Get the corresponding exon ID

            num_samples = mus_pids[pid]['num-samples'][target_index];

            target_tracker[target] = { 'gid' : mus_pids[pid]['gid'], 'tid' : mus_pids[pid]['tid'], 'eid' : meid, 'pid' : pid,
                                        'init sample hits' : str(num_samples), 'exonerate sample hits' : 0 };
            #target_tracker[target] = { 'gid' : mus_pids[pid]['gid'], 'tid' : mus_pids[pid]['tid'], 'eid' : meid, 'pid' : pid,
            #                            'exonerate sample hits' : 0 };
            # Initialize the target tracker for the current target.

            if debug:
                print(target, ":", num_samples, filter_thresh);

            if num_samples < filter_thresh:
                continue;
            # THE FILTER
            
            exonerate_file = os.path.join(indir, target + "-exonerate.txt");

            if not os.path.isfile(exonerate_file):
                continue;
            # Retrieve the current exonerate output file.

            sequence_tracker['mm10']['nt-seq'] += mus_pids[pid]['mm-nt-seqs'][target_index];
            sequence_tracker['mm10']['aa-seq'] += mus_pids[pid]['mm-aa-seqs'][target_index];

            if debug:
                print(target, ":", exonerate_file);

            sample_segment_tracker[target] = {};
            # Add the target to the main tracker dict.

            gff_flag, seq_flag, mouse_seq_flag, aln_flag = False, False, False, False;
            # The exonerate output files have a lot o information in chunks, including:
            # 1. GFF info for each alignment.
            # 2. Aligned CDS from the target.
            # 3. Aligned mouse amino acid sequence.
            # 4. Info about the alignment.
            # 5. The alignment (1 site per line).
            # These flags help track which chunk we're reading through.

            ##########
            # Begin line parsing
            for line in open(exonerate_file):
                if any(line.startswith(skip_str) for skip_str in ["Command line", "Hostname", "#", "--"]):
                    continue;
                # Skip lines that start with some key characters
                # if pid == "ENSMUSP00000021056" and cur_sample == "Tarsomys_apoensis_KU167148":
                #     print(line);
                #     print(gff_flag, seq_flag, mouse_seq_flag, aln_flag);

                ##########
                # Begin GFF block parse
                if line == "# --- START OF GFF DUMP ---\n":
                    gff_flag = True;
                    continue;

                if gff_flag:
                    if line == "# --- END OF GFF DUMP ---\n":
                        gff_flag = False;
                    continue;
                # For now, none of the GFF info is used, but its good to have in the output files.
                # End GFF block parse
                ##########

                ##########
                # Begin INFO block parse
                # The INFO line contains information about the current alignment and is the basis for adding targets.
                # Since it is only one line there's no need for a flag.

                split_line = line.strip().split("\t");
                # Split the line by tab to access info.
                if split_line[0] == "INFO":
                    cur_sample = split_line[1];
                    # Get the name of the target sample.

                    # if pid == "ENSMUSP00000021056" and cur_sample == "Tarsomys_apoensis_KU167148":
                    #     print(split_line);

                    target_strand = split_line[2];
                    if target_strand not in ["+", "-"]:
                        print("INVALID TARGET STRAND: " + target_strand);
                        print(mtid, pid, meid, target, cur_sample, line);
                        return "exit";                             

                    if cur_sample not in sample_segment_tracker[target]:
                        sample_segment_tracker[target][cur_sample] = { 'gid' : mus_pids[pid]['gid'], 'tid' : mus_pids[pid]['tid'], 'eid' : meid, 'pid' : pid,
                            'seg-start' : [], 'seg-end' : [], 'seg-len' : [], 'seg-sim' : [], 'seg-div3' : [], 'seg-stops' : [], 
                            'num-introns-pred' : [], 'num-segs-pred' : 0,
                            'nt-seq' : [], 'aa-seq' : [], 'mm-nt-seq' : [], 'mm-aa-seq' : [],
                            'nt-aln' : [], 'aa-aln' : [], 'mm-nt-aln' : [], 'mm-aa-aln' : [] };

                    if cur_sample not in sample_target_tracker[pid]: 
                        sample_target_tracker[pid][cur_sample] = { 'gid' : mus_pids[pid]['gid'], 'tid' : mus_pids[pid]['tid'], 'eid' : meid, 'pid' : pid,
                                                                    'mm-exons' : len(mus_pids[pid]['eids']), 'exonerate hits' : 0 }

                    if cur_sample not in sequence_tracker:
                        sequence_tracker[cur_sample] = {'nt-seq' : "", 'aa-seq' : ""};

                    if split_line[2] == "+":
                        start, end, length, sim = int(split_line[3]), int(split_line[4]), int(split_line[5]), float(split_line[8]);
                    elif split_line[2] == "-":
                        start, end, length, sim = int(split_line[4]), int(split_line[3]), int(split_line[5]), float(split_line[8]);
                    # Retrieve the info for the current alignment. The only difference between strands is the order of the start
                    # and end coordinates. Exonerate handles everything else.

                    sample_segment_tracker[target][cur_sample]['seg-start'].append(start);
                    sample_segment_tracker[target][cur_sample]['seg-end'].append(end);
                    sample_segment_tracker[target][cur_sample]['seg-len'].append(length);
                    sample_segment_tracker[target][cur_sample]['seg-sim'].append(sim);                   
                    sample_segment_tracker[target][cur_sample]['num-introns-pred'].append(0);
                    sample_segment_tracker[target][cur_sample]['num-segs-pred'] += 1;
                    # Add info to the sample_segment_tracker dict.

                    #print(cur_sample);
                # End INFO block parse
                ##########

                ##########
                # Begin target sequence block parse
                if line == "BEGIN SEQUENCE\n":
                    seq_flag = True;
                    cur_seq = "";
                    continue;
                # Start the sequence block and initialize the sequence string.

                if seq_flag:
                    if line == "END SEQUENCE\n":
                        sample_segment_tracker[target][cur_sample]['nt-seq'].append(cur_seq);
                        sample_segment_tracker[target][cur_sample]['aa-seq'].append(str(Seq(cur_seq).translate()));

                        seq_flag = False;
                        continue;
                    # Once the END SEQUENCE string is encountered, save the built up nucleotide sequence and the translation.

                    cur_seq += line.strip();
                    # Add on the current line of sequence to the cur_seq string.
                # This block retrieves the aligned coding sequence as determined by exonerate. This is currently
                # unused in lieu of manual parsing of the alignment sequence below.
                # Begin target sequence block parse
                ##########
                
                ##########
                # Begin mouse sequence block parse
                if line == "BEGIN MOUSE\n":
                    mouse_seq_flag = True;
                    cur_seq = "";
                    continue;
                # Start the sequence block and initialize the sequence string.

                if mouse_seq_flag:
                    if line == "END MOUSE\n":
                        sample_segment_tracker[target][cur_sample]['mm-aa-seq'].append(cur_seq);
                        #target_tracker[pid][target][cur_sample]['mm-aa-seq'].append(str(Seq(cur_seq).translate()));
                        mouse_seq_flag = False;
                        continue;
                    # Once the END MOUSE string is encountered, save the built up aa sequence.

                    cur_seq += line.strip();
                    # Add on the current line of sequence to the cur_seq string.
                # This block retrieves the mouse protein sequence that was aligned to the target. This is used later on
                # to determine similarity of chunks of target sequence around stop codons.
                # Begin mouse sequence block parse
                ##########

                ##########
                # Begin alignment block parse
                if line.startswith("BEGIN ALN"):
                    aln_flag = True;
                    cur_mm_aln = "";
                    cur_aln = "";
                    continue;
                # Start the alignment block and initialize the sequence strings.

                if aln_flag:
                    if line == "END ALN\n":
                        # print(cur_sample);
                        # print(cur_aln);
                        # print(cur_mm_aln);
                        sample_segment_tracker[target][cur_sample]['nt-aln'].append(cur_aln);
                        sample_segment_tracker[target][cur_sample]['aa-aln'].append(str(Seq(cur_aln).translate()));
                        sample_segment_tracker[target][cur_sample]['mm-aa-aln'].append(cur_mm_aln);
                        #target_tracker[pid][target][cur_sample]['mm-aa-aln'].append(str(Seq(cur_mm_aln).translate()));
                        aln_flag = False;
                        continue;
                    # Once the END ALN string is encountered, save the built up nucleotide alignment and its translation and the corresponding mouse amino acid sequence.

                    #print(line);
                    line = line.strip().split("\t");
                    aa, aa_adv, codon, codon_adv, aln_type, aln_label = line[0], int(line[1]), line[4], int(line[5]), line[8], line[10];
                    # Split the line and save the relevant information.
                    
                    if aln_type not in aln_type_str:
                        print("UNKONWN ALIGN TYPE: " + aln_type);
                        print(mtid, pid, meid, target, cur_sample, line);
                        return "exit";
                    # Check if the alignment type is one of the known exonerate types. If not, report and exit.

                    if aln_type in ["start to match", "match to end"]:
                        continue;
                    # The start and end of the alignments are denoted by gaps in both sequence with these labels. Skip these.

                    if aln_label in ["5'ss", "intron", "3'ss", "none"]:
                        continue;
                    # Several types of lines don't require anything, including: splice sites at introns, intron sites, and gaps in both sequences (none). Skip.

                    if aln_type in intron_start_str:
                        target_tracker[pid][target][cur_sample]['num-introns'][-1] += 1;
                        continue;
                    # If alignment type is intron start, add to count of introns for target segment and skip.

                    if aln_type in intron_split_str:
                        split_codon_part = codon;
                        continue;
                    # If an intron starts with a split codon, keep track of the codon part and skip. All other intron lines should be skipped until
                    # end based on above.

                    if aln_type in intron_end_str:
                        codon = split_codon_part + codon;
                        if debug:
                            print("INTRON END:", mtid, pid, meid, target, cur_sample, line, split_codon_part, codon);
                        split_codon_part = "";
                    # At the end of an intron, if the codon was split, add the parts together here.

                    if aa_adv not in [0, 1]:                        
                        print("Unknown mouse AA advance:");
                        print(pid, target, cur_sample, line);
                        return "exit";
                    # # The mouse sequence should only ever jump ahead by these intervals. If not, report and exit.

                    if codon_adv not in [0, 1, 2, 3]:
                        print("Unknown target codon advance:");
                        print(pid, target, cur_sample, line);
                        return "exit";
                    # # The target sequence should only ever jump ahead by these intervals. If not, report and exit.

                    if len(codon) < 3:
                        codon = codonComplete(codon);
                        if debug:
                            print("FRAMESHIFT/INSERTION:", mtid, pid, meid, target, cur_sample, line, codon);
                    # For gaps in the target sequence or frameshifts, exonerate only reports partial codons. These must be expanded for the final sequence.

                    cur_aln += codon;
                    cur_mm_aln += aa;
                    # Add the current sequences to the sequence strings.
                # End alignment block parse
                ##########
            # End line parsing
            ##########     
      
            # if target == "CDS_chr11_id00424":
            #     return "exit";
            
            ##########
            # Begin segment selection block
            for s in sample_segment_tracker[target]:
            # For every sample in the current target, we need to select the longest predicted sequence, deal with stop codons, and
            # concatenate the sequences together in the sample_tracker.

                num_segs = sample_segment_tracker[target][s]['num-segs-pred'];
                # Get the number of predicted segments for the current target/sample combo.

                if num_segs < 1:
                    continue;
                # If no segments are predicted, skip. This shouldn't happen?

                for seg_ind in range(len(sample_segment_tracker[target][s]['seg-start'])):
                    cur_nt_len = sample_segment_tracker[target][s]['seg-len'][seg_ind];
                    cur_aa_aln = sample_segment_tracker[target][s]['aa-aln'][seg_ind];

                    if cur_nt_len % 3 == 0:
                        sample_segment_tracker[target][s]['seg-div3'].append("TRUE");
                    else:
                        sample_segment_tracker[target][s]['seg-div3'].append("FALSE");
                    sample_segment_tracker[target][s]['seg-stops'].append(cur_aa_aln.count("*"))
                # Get some info about frame and stop codons for each predicted segment.

                if debug:
                    print(mtid, pid, meid, target, s);

                max_len_ind = 0;
                if num_segs > 1:
                    max_len = max(sample_segment_tracker[target][s]['seg-len']);
                    max_len_ind = sample_segment_tracker[target][s]['seg-len'].index(max_len); 
                # If more than one segment is predicted for this target/sample combo, get the index of the longest segment. If
                # only one segment is predicted, the index is set at 0 to retrieve the only prediction in each list (at index 0).

                sample_segment_tracker[target][s]['seg-start'] = sample_segment_tracker[target][s]['seg-start'][max_len_ind];
                sample_segment_tracker[target][s]['seg-end'] = sample_segment_tracker[target][s]['seg-end'][max_len_ind];
                sample_segment_tracker[target][s]['seg-len'] = sample_segment_tracker[target][s]['seg-len'][max_len_ind];
                sample_segment_tracker[target][s]['seg-sim'] = sample_segment_tracker[target][s]['seg-sim'][max_len_ind];
                sample_segment_tracker[target][s]['seg-div3'] = sample_segment_tracker[target][s]['seg-div3'][max_len_ind];
                sample_segment_tracker[target][s]['num-introns-pred'] = sample_segment_tracker[target][s]['num-introns-pred'][max_len_ind];
                sample_segment_tracker[target][s]['seg-stops'] = sample_segment_tracker[target][s]['seg-stops'][max_len_ind];
                # Save info from longest predicted segments.

                sample_segment_tracker[target][s]['nt-seq'] = sample_segment_tracker[target][s]['nt-seq'][max_len_ind];
                sample_segment_tracker[target][s]['aa-seq'] = sample_segment_tracker[target][s]['aa-seq'][max_len_ind];
                sample_segment_tracker[target][s]['mm-aa-seq'] = sample_segment_tracker[target][s]['mm-aa-seq'][max_len_ind];
                # Save the longest predicted sequences.

                cur_nt_aln = sample_segment_tracker[target][s]['nt-aln'][max_len_ind];
                sample_segment_tracker[target][s]['nt-aln'] = cur_nt_aln;
                cur_aa_aln = sample_segment_tracker[target][s]['aa-aln'][max_len_ind];
                sample_segment_tracker[target][s]['aa-aln'] = cur_aa_aln;
                cur_mm_aa_aln = sample_segment_tracker[target][s]['mm-aa-aln'][max_len_ind];
                sample_segment_tracker[target][s]['mm-aa-aln'] = cur_mm_aa_aln;
                # Save the corresponding alignments from the longest predicted sequences.

                #target_tracker[pid][target][s]['orig-seg-sim'] = percIdent(cur_aa_aln, cur_mm_aa_aln) * 100;
                # Save the alignment similarity of the longest predicted sequence.

                num_stops = sample_segment_tracker[target][s]['seg-stops'];
                # Count the number of stop codons in the longest predicted alignment and save.

                ##########
                # Start stop codon block
                do_stop = False;
                if num_stops > 0 and do_stop:
                    # If there are stop codons in the longest predicted sequence, we need to split the alignment and select the one with
                    # the highest similarity to mouse.

                    aa_chunks = cur_aa_aln.split("*");
                    # Split the amino acid alignment by stop codon symbol.

                    stop_inds = [ i for i, residue in enumerate(cur_aa_aln) if residue == "*" ];
                    # Get a list of indices of stop codons in the amino acid alignment.

                    nt_stop_inds = [ i*3 for i in stop_inds ];
                    # Adjust the indices for the nucleotide alignment.

                    stop_codons = [];
                    for ind in nt_stop_inds:
                        stop_codons.append(cur_nt_aln[ind:ind+3]);
                    # Get the stop codons to check in a sec.

                    mm_aa_chunks = splitSeq(cur_mm_aa_aln, stop_inds, "aa");
                    # Split the mouse amino acid alignment by the stop codon indices.

                    nt_chunks = splitSeq(cur_nt_aln, nt_stop_inds, "codon");
                    # Split the target nucleotide alignment by the stop codon indices.

                    mm_nt_chunks = splitSeq(cur_mm_nt_aln, nt_stop_inds, "codon");

                    # print(aa_chunks);
                    # print(mm_chunks);
                    # print(nt_chunks);
                    # print('---')

                    if "**" in cur_aa_aln or cur_aa_aln[0] == "*" or cur_aa_aln[-1] == "*":
                        aa_chunks = list(filter(None, aa_chunks));
                        mm_aa_chunks = list(filter(None, mm_aa_chunks));
                        nt_chunks = list(filter(None, nt_chunks));
                        mm_nt_chunks = list(filter(None, mm_nt_chunks));
                    # Sequences with consecutive stop codons leave an empty string in the chunk lists. Remove them here.

                    # print(aa_chunks);
                    # print(mm_chunks);
                    # print(nt_chunks);

                    target_tracker[pid][target][s]['stop-chunk-lens'] = [ len(chunk) for chunk in aa_chunks ];

                    min_chunk_len = 5;
                    if all(clen < 5 for clen in target_tracker[pid][target][s]['stop-chunk-lens']):
                        min_chunk_len = 0;

                    if debug:
                        print("------");
                        print("STOP CODON REPORT");
                        print("\t", mtid, pid, meid, target, s);
                        print("\tOriginal NT aln      : ", cur_nt_aln);
                        print("\tStop codons         : ", stop_codons);
                        print("\tSplit NT aln         : ", nt_chunks);
                        print("\tOriginal mouse NT aln: ", cur_mm_nt_aln);
                        print("\tSplit mouse NT aln   : ", mm_nt_chunks);
                        print("\t----");
                        print("\tOriginal AA aln      : ", cur_aa_aln);
                        print("\tSplit AA aln         : ", aa_chunks);
                        print("\tOriginal mouse AA aln: ", cur_mm_aa_aln);
                        print("\tSplit mouse AA aln   : ", mm_aa_chunks);
                        print("\tAA stop codon indices: ", stop_inds);
                        print("\tNT stop codon indices: ", nt_stop_inds);
                        print("\t----");

                    chunk_sims = [ percIdent(aa_chunks[i], mm_aa_chunks[i]) * 100  for i in range(len(aa_chunks)) ];
                    #for i in range(len(aa_chunks)):
                    #    chunk_sims.append(SequenceMatcher(None, aa_chunks[i], mm_chunks[i],).ratio() * 100);
                    # For each chunk of sequence around a stop codon, calculate amino acid sequence similarity between mouse and target.

                    if debug:
                        print("\tMin chunk len        : ", min_chunk_len);
                        print("\tLength of chunks     : ", target_tracker[pid][target][s]['stop-chunk-lens']);
                        print("\tOriginal similarity  : ", target_tracker[pid][target][s]['orig-seg-sim']);
                        print("\tSimilarity of chunks : ", chunk_sims);
                        print("\t------");

                    assert all(c in ["TAG", "TAA", "TAR", "TGA", "TRA"] for c in stop_codons), "\nINVALID STOP CODON FOUND";

                    target_tracker[pid][target][s]['stop-chunk-sims'] = ",".join( [ str(s) for s in chunk_sims ] );
                    # Save the list of similarities for each chunk.

                    chunk_sims = [ chunk_sims[i] if len(aa_chunks[i]) >= min_chunk_len else 0.0 for i in range(len(chunk_sims)) ];

                    max_sim = max(chunk_sims);
                    max_sim_ind = chunk_sims.index(max_sim);
                    target_tracker[pid][target][s]['final-seg-sim'] = str(max_sim);
                    # Select the chunk with the highest similiarity and save.

                    if debug:
                        print("\tMax similarity       : ", max_sim);
                        print("\t------");

                    target_tracker[pid][target][s]['nt-aln'] = nt_chunks[max_sim_ind];
                    target_tracker[pid][target][s]['aa-aln'] = aa_chunks[max_sim_ind];
                    # Get the alignment segments with the max similarity

                    if debug:
                        print("\tFinal AA aln         : ", target_tracker[pid][target][s]['aa-aln']);
                        print("\tFinal mouse AA aln   : ", mm_aa_chunks[max_sim_ind]);
                        print("\tFinal NT aln         : ", target_tracker[pid][target][s]['nt-aln']);
                        print("\tFinal mouse NT aln   : ", mm_nt_chunks[max_sim_ind]);
                        print("\t------");

                    # if pid == "ENSMUSP00000091532" and target == "CDS_chr11_id00019" and s == "Sundamys_muelleri_JAE1283":
                    #     return "exit";
                    #return "exit";
                # End stop codon block
                ##########

                sample_segment_tracker[target][s]['nt-seq'] = sample_segment_tracker[target][s]['nt-aln'].replace("-", "");
                sample_segment_tracker[target][s]['aa-seq'] = sample_segment_tracker[target][s]['aa-aln'].replace("-", "");
                # Save the UNALIGNED sequences now by replacing gaps with empty strings in the alignments.

                sample_segment_tracker[target][s]['seg-len'] = len(sample_segment_tracker[target][s]['nt-seq']);
                # Save the length of the final sequence.

                sample_segment_tracker[target][s]['seg-sim'] = percIdent(sample_segment_tracker[target][s]['aa-aln'], sample_segment_tracker[target][s]['mm-aa-aln']) * 100;
                # Recalculate the alignment similarity.

                # target_tracker[pid][target][s]['final-seg-div3'] = "TRUE";
                # if target_tracker[pid][target][s]['final-seg-len'] % 3 != 0:
                #     target_tracker[pid][target][s]['final-seg-div3'] = "FALSE";
                # # Determine if the final segment is divisible by 3.

                sequence_tracker[s]['nt-seq'] += sample_segment_tracker[target][s]['nt-seq'];
                sequence_tracker[s]['aa-seq'] += sample_segment_tracker[target][s]['aa-seq'];
                # Finally, concatenate the final predicted sequence to the sequence in the sample_tracker for output later.

                sample_target_tracker[pid][s]['exonerate hits'] += 1;
                target_tracker[target]['exonerate sample hits'] += 1;
            # End segment selection block
            ##########
            # if target == "CDS_chr11_id00010":
            #     return target_tracker, written, skipped;
            
        # End target block
        ##########   


        ##########
        # Begin output sequence block
        if noseq:
            skipped += 1;
            continue;
        else:       
            nt_outfile = os.path.join(outdir, "nt", pid + ".fa");
            aa_outfile = os.path.join(outdir, "aa", pid + ".fa");
            # Define the output file name for the current protein ID

            with open(nt_outfile, "w") as ntout, open(aa_outfile, "w") as aaout:
            # Open both the CDS and AA output files.
                #ntout.write(">mm10\n" + mus_pids[pid][] + "\n");
                #aaout.write(">mm10\n" + str(Seq(mseq[mtid]).translate()).replace("*", "") + "\n");
                # Retrieve the mouse CDS sequence, translate, and write to files.

                # add_rat = False;
                # if rtid in rseq:
                #     add_rat = True;
                #     ntout.write(">rnor6\n" + rseq[rtid] + "\n");
                #     aaout.write(">rnor6\n" + str(Seq(rseq[rtid]).translate()).replace("*", "") + "\n");
                # Determine if there is an orthologous rat sequence to add and add here.

                ##########
                # Begin sample output block
                for sample in sequence_tracker:
                # Write sequences for every sample that has a sequence predicted for this target.
                    ntout.write(">" + sample + "\n" + sequence_tracker[sample]['nt-seq'] + "\n");
                    aaout.write(">" + sample + "\n" + sequence_tracker[sample]['aa-seq'] + "\n");
                    # Write both the sequences.
    
                    # if sample == "mm10":
                    #     continue;

                    # cur_div3 = "TRUE";
                    # cur_tot_len = len(sample_tracker[sample]['nt-seq']);
                    # if cur_tot_len % 3 != 0:
                    #     cur_div3 = "FALSE";
                    # # Determine if the current nucleotide sequence is divisble by 3.

                    # cur_stops = sample_tracker[sample]['aa-seq'].count("*");
                    # # Determine if there are any stop codons in the final sequence.                

                    # for target in sample_tracker[sample]['targets']:
                    #     target_tracker[pid][target][sample]['tot-cds-len'] = cur_tot_len;
                    #     target_tracker[pid][target][sample]['cds-div3'] = cur_div3;
                    #     target_tracker[pid][target][sample]['cds-stops'] = cur_stops;
                    #     target_tracker[pid][target][sample]['num-cds-protein'] = len(sample_tracker[sample]['targets']);
                    # Add te final CDS info for every target in the sample for the current protein.
                # End sample output block
                ##########

            #target_tracker[pid]["NA"] = {};
            #target_tracker[pid]["NA"]["mm10"] = { 'num-cds-protein' : str(mexp[mtid]['num-cds']), 'num-segs-target' : "NA", 'num-introns' : "NA", 'seg-start' : "NA", 'seg-end' : "NA", 
            #                                'tot-cds-len' : str(mexp[mtid]['tot-len']), 'orig-seg-len' : "NA", 'cds-div3' : "NA", 'orig-seg-div3' : "NA", 'orig-seg-sim' : "NA", 
            #                                'cds-stops' : "NA", 'seg-stops' : "NA",  'nt-seq' : "NA", 'aa-seq' : "NA", 'mm-seq' : "NA", 'nt-aln' : "NA", 'aa-aln' : "NA", 'mm-aln' : "NA", 
            #                                'stop-chunk-sims' : "NA", 'stop-chunk-lens' : "NA", 'final-seg-len' : "NA", 'final-seg-sim' : "NA", "final-seg-div3" : "NA" };
            # if add_rat:
            #     target_tracker[pid]["NA"]["rnor6"] = { 'num-cds-protein' : str(rexp[rtid]['num-cds']), 'num-segs-target' : "NA", 'num-introns' : "NA", 'seg-start' : "NA", 'seg-end' : "NA", 
            #                                 'tot-cds-len' : str(rexp[rtid]['tot-len']), 'orig-seg-len' : "NA", 'cds-div3' : "NA", 'orig-seg-div3' : "NA", 'orig-seg-sim' : "NA", 
            #                                 'cds-stops' : "NA", 'seg-stops' : "NA",  'nt-seq' : "NA", 'aa-seq' : "NA", 'mm-seq' : "NA", 'nt-aln' : "NA", 'aa-aln' : "NA", 'mm-aln' : "NA", 
            #                                 'stop-chunk-sims' : "NA", 'stop-chunk-lens' : "NA", 'final-seg-len' : "NA", 'final-seg-sim' : "NA", "final-seg-div3" : "NA" };
            # Add the info for mouse and rat by protein.
                
            written += 1;
            # End output sequence block
            ##########
    # End protein block
    ##########

    return target_tracker, sample_target_tracker, sample_segment_tracker, written, skipped;

###########################################################
# Options

parser = argparse.ArgumentParser(description="Exonerate command generator");
parser.add_argument("-i", dest="input", help="Directory of input exonerate outputs.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for un-aligned, concatenated CDS seqs.", default=False);
parser.add_argument("-f", dest="filter", help="Filter threshold. Exons with fewer than this many matching samples won't be retrieved.", type=int, default=0);
parser.add_argument("-p", dest="procs", help="Number of processes to use. Default: 1", type=int, default=1);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--noseq", dest="noseq", help="Set to not write the sequences, just make the counts.", action="store_true", default=False);
args = parser.parse_args();
# IO options

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory must be defined with -i.");
args.input = os.path.abspath(args.input);

if not args.output:
    sys.exit( " * Error 2: An output directory must be defined with -o.");

if args.filter < 0:
    sys.exit( " * Error 3: Filter threshold (-f) must be a positive integer or zero.");

if args.procs < 1:
    sys.exit( " * Error 4: Number of processes (-p) must be a positive integer.");

args.output = os.path.abspath(args.output);
if args.output[-1] == "/":
    args.output = args.output[:-1];
args.output += "-f" + str(args.filter);
# Adjust the output directory name to add the filter threshols.

if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 5: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

info_pad = 26
pad = 55;
cwd = os.getcwd();
# Job vars

annotation_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/repro_exons.selected_txpts.exon_counts.txt"
# Table matching protein, transcript, gene, and exon IDs

mouse_nt_dir = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/repro_genes-selected-cds-nt-trimmed/"
mouse_aa_dir = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/repro_genes-selected-cds-aa-trimmed/"
# The mouse trimmed exon sequences.

if not os.path.isdir(args.output):
    print("# Creating output directory.");
    os.system("mkdir " + args.output);
# Create the output directory.

log_file = os.path.join(args.output, "exonerate-to-cds.log");
segments_file = os.path.join(args.output, "sample-segments-exonerate.csv");
samples_file = os.path.join(args.output, "sample-targets-exonerate.csv");
targets_file = os.path.join(args.output, "targets-exonerate.csv");
# The various output files.
    
#"logs", os.path.basename(os.path.normpath(args.output)) + "-exonerate-to-cds-indiv.log");
# Job files

##########################
# Reporting run-time info for records.

with open(log_file, "w") as logfile:
    core.runTime("#!/bin/bash\n# Exonerate sequence extractor", logfile);
    core.PWS("# IO OPTIONS", logfile);
    core.PWS(core.spacedOut("# Input directory:", info_pad) + args.input, logfile);
    core.PWS(core.spacedOut("# Annotation file:", info_pad) + annotation_file, logfile);
    core.PWS(core.spacedOut("# Mouse NT seqs:", info_pad) + mouse_nt_dir, logfile);
    core.PWS(core.spacedOut("# Mouse AA seqs:", info_pad) + mouse_aa_dir, logfile);
    core.PWS(core.spacedOut("# Sample filter threshold:", info_pad) + str(args.filter), logfile);
    core.PWS(core.spacedOut("# Output directory:", info_pad) + args.output, logfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", info_pad) + "Overwriting previous files in output directory.", logfile);
    if args.noseq:
        core.PWS("# --noseq set. Not writing sequences.", logfile);
    else:
        os.system("mkdir " + args.output + "/nt/");
        os.system("mkdir " + args.output + "/aa/");
    core.PWS(core.spacedOut("# Segments file:", info_pad) + segments_file, logfile);
    core.PWS(core.spacedOut("# Samples file:", info_pad) + samples_file, logfile);
    core.PWS(core.spacedOut("# Targets file:", info_pad) + targets_file, logfile);
    core.PWS(core.spacedOut("# Procs:", info_pad) + str(args.procs), logfile);
    core.PWS(core.spacedOut("# Logfile:", info_pad) + log_file, logfile);
    core.PWS("# ----------------", logfile);
    # Log info

    ##########################
    # Reading annotation info.
    #transcripts_file = "../02-Annotation-data/selected-transcripts-targets.txt";
    # The selected transcripts file based on the mouse and rat references.

    #mus_transcripts = "../02-Annotation-data/mm10-cds.tab";
    #mus_seq_file = "../Reference-genomes/mm10/Mus_musculus.GRCm38.cds.all.fa";
    # Mouse sequences and annotation info. The tab file is from gxf_parse.py.

    #rat_transcripts = "../02-Annotation-data/rnor6-cds.tab";
    #rat_seq_file = "../Reference-genomes/Rnor6/Rattus_norvegicus.Rnor_6.0.cds.all.fa";
    #rat_id_file = "../Reference-genomes/Rnor6/rnor6-ens99-ids.tab";
    # Rat sequences, annotation info, and ID conversions. The tab file is from gxf_parse.py.
    # Various annotation files...

    skip_file = "mm10-exons-w-stops.txt";
    #skip_eids = [ line.strip().split(",")[1] for line in open(skip_file, "r") ];
    skip_ids = [ line.strip() for line in open(skip_file, "r") ];

    core.PWS(core.spacedOut("# " + core.getDateTime() + " Reading main annotation file:", pad) + annotation_file, logfile);
    mus_pids = {};
    # The mus_pids dict contains the conversion between mouse IDs and target/contig IDs.
    num_targets = 0;
    first = True;
    for line in open(annotation_file):
        if first:
            first = False;
            continue;
        #line = line.replace("\"", "").strip().split(",");
        line = line.replace("\"", "").strip().split("\t");
        #targ, gid, tid, pid, eid, num_samples = line[0].replace("\"", ""), line[5].replace("\"", ""), line[9].replace("\"", ""), line[13].replace("\"", ""), line[17].replace("\"", ""), line[22];
        targ, egn, gid, tid, pid, eid, num_samples = line[0].replace("\"", "")+"_"+line[4].replace("\"", ""), line[0].replace("\"", ""), line[1].replace("\"", ""), line[2].replace("\"", ""), line[3].replace("\"", ""), line[4].replace("\"", ""), line[13];
        pid_eid = ",".join((pid, eid));
        if pid == "NA" or pid_eid in skip_ids:
            continue;

        if pid not in mus_pids:
            mus_pids[pid] = { 'gid' : gid, 'tid' : tid, 'eids' : [], 'targets' : [], 'num-samples' : [], 'mm-nt-seqs' : [], 'mm-aa-seqs' : [] };
            #mus_pids[pid] = { 'gid' : gid, 'tid' : tid, 'eids' : [], 'targets' : [], 'mm-nt-seqs' : [], 'mm-aa-seqs' : [] };

        mus_pids[pid]['eids'].append(eid);
        mus_pids[pid]['targets'].append(targ);
        mus_pids[pid]['num-samples'].append((int(num_samples)));

        #nt_file = glob(os.path.join(mouse_nt_dir, eid + "chr*.fa"))[0];
        nt_file = os.path.join(mouse_nt_dir, eid + "-nt-trimmed.fa");
        cur_title, cur_seq = mseq.fastaGetLists(nt_file);
        mus_pids[pid]['mm-nt-seqs'].append(cur_seq[0]);
        
        #aa_file = glob(os.path.join(mouse_aa_dir, "translated_" + eid + "chr*.fa-out.fas-out.fas"))[0];
        aa_file = os.path.join(mouse_aa_dir, eid + "-aa-trimmed.fa");
        cur_title, cur_seq = mseq.fastaGetLists(aa_file);
        mus_pids[pid]['mm-aa-seqs'].append(cur_seq[0]);

        num_targets += 1;
    core.PWS(core.spacedOut("# Protein IDs read:", info_pad) + str(len(mus_pids)), logfile);
    core.PWS(core.spacedOut("# Total targets:", info_pad) + str(num_targets));
    core.PWS("# ----------------", logfile);
    # Reads the mouse and target annotation info.

      
    # core.PWS("# Reading mouse seqs: " + mus_seq_file + " " + mus_transcripts, logfile);
    # mus_seqs, mus_expected = readSeqInfo(mus_seq_file, mus_transcripts);
    # core.PWS("# ----------------", logfile);
    # # Reads and stores the mouse sequences.

    # core.PWS("# Reading rat seqs: " + rat_seq_file + " " + rat_transcripts, logfile);
    # rat_seqs, rat_expected = readSeqInfo(rat_seq_file, rat_transcripts);
    # core.PWS("# ----------------", logfile);
    # # Reads and stores the rat sequences.
    
    # core.PWS("# Reading selected transcripts: " + transcripts_file, logfile);
    # mus_pid_to_rat_pid = {};
    # first = True;
    # for line in open(transcripts_file):
    #     if line[0] == "#":
    #         continue;
    #     if first:
    #         first = False;
    #         continue;
    #     # Skip the comment lines and first line.

    #     line = line.strip().split("\t");
    #     mus_pid, rat_pid = line[2], line[6];
    #     # Split and unpack the line info.

    #     mus_pid_to_rat_pid[mus_pid] = rat_pid;
    #     # Add to conversion list.
    # core.PWS("# IDs read: " + str(len(mus_pid_to_rat_pid)), logfile);
    # core.PWS("# ----------------", logfile);
    # # Read the selected transcripts file to get converions between mouse and rat protein IDs.

    # core.PWS("# Reading rat IDs: " + rat_id_file, logfile);
    # rat_pid_to_tid = {};
    # for line in open(rat_id_file):
    #     line = line.strip().split("\t");
    #     if len(line) != 3:
    #         continue;
    #     rat_pid_to_tid[line[2]] = line[1];
    # core.PWS("# Rat IDs read: " + str(len(rat_pid_to_tid)), logfile);
    # core.PWS("# ----------------", logfile);
    # # Read the rat ID file to get conversions between rat protein and transcript IDs.

    ##########################
    # Main parsing block -- splits out into multiple processes.

    written_main, skipped_main = 0,0;
    tracker_main = {};
    target_tracker_main = {};
    sample_target_tracker_main = {};
    sample_segment_tracker_main = {};
    # The main count and tracker variables. The trackers will be updated for each chunk of proteins.

    pids = list(mus_pids.keys());
    pids_per_proc = int(math.ceil(len(pids) / args.procs));
    pid_lists = list(core.chunks(pids, pids_per_proc));
    # Splitting the protein IDs into even chunks for parallelization

    #pool = mp.Pool(processes=args.procs);
    # Starting the pool of processes

    core.PWS(core.spacedOut("# " + core.getDateTime() + " Start exonerate parse.", pad), logfile);
    #for result in pool.starmap(parseExonerate, ( (pl, mus_pids, args.input, args.output, args.filter, args.noseq, True) for pl in pid_lists )):

    for pl in pid_lists:
        result = parseExonerate(pl, mus_pids, args.input, args.output, args.filter, args.noseq, True);
        if result == "exit":
            sys.exit();
        target_tracker_main.update(result[0]);
        sample_target_tracker_main.update(result[1]);
        sample_segment_tracker_main.update(result[2]);
        written_main += result[3];
        skipped_main += result[4];
    # Pass each chunk of IDs to the parse function and collect variables
    core.PWS(core.spacedOut("# " + core.getDateTime() + " End exonerate parse.", pad), logfile);
    core.PWS(core.spacedOut("# Protein seqs written:", info_pad) + str(written_main), logfile);
    core.PWS(core.spacedOut("# Protein seqs skipped:", info_pad) + str(skipped_main), logfile);
    core.PWS("# ----------------", logfile);

    ##########################
    # Output.
    core.PWS(core.spacedOut("# " + core.getDateTime() + " Writing segments file:", pad) + segments_file, logfile);
    with open(segments_file, "w") as outfile:
        first = True;
        for target in sample_segment_tracker_main:
            for sample in sample_segment_tracker_main[target]:
                if first:
                    givens = ["target", "sample"];
                    headers = givens + list(sample_segment_tracker_main[target][sample].keys());
                    skippers = ["aa-seq","nt-seq","nt-aln","aa-aln",'mm-aa-seq','mm-nt-seq','mm-nt-aln','mm-aa-aln'];
                    headers = [h for h in headers if h not in skippers];
                    outfile.write(",".join(headers) + "\n");
                    first = False;

                outline = [target, sample] + [ str(sample_segment_tracker_main[target][sample][col]) for col in headers if col not in givens ];
                outfile.write(",".join(outline) + "\n");

    core.PWS(core.spacedOut("# " + core.getDateTime() + " Writing samples file:", pad) + samples_file, logfile);
    with open(samples_file, "w") as outfile:
        first = True;
        for pid in sample_target_tracker_main:
            for sample in sample_target_tracker_main[pid]:
                if first:
                    givens = ["pid", "sample"];
                    headers = givens + list(sample_target_tracker_main[pid][sample].keys());
                    outfile.write(",".join(headers) + "\n");
                    first = False;

                outline = [pid, sample] + [ str(sample_target_tracker_main[pid][sample][col]) for col in headers if col not in givens ];
                outfile.write(",".join(outline) + "\n");

    core.PWS(core.spacedOut("# " + core.getDateTime() + " Writing targets file:", pad) + targets_file, logfile);
    with open(targets_file, "w") as outfile:
        first = True;
        for target in target_tracker_main:
            if first:
                givens = ["target"];
                headers = givens + list(target_tracker_main[target].keys());
                outfile.write(",".join(headers) + "\n");
                first = False;

            outline = [target] + [ str(target_tracker_main[target][col]) for col in headers if col not in givens ];
            outfile.write(",".join(outline) + "\n");   



    # headers = ["pid","sample","target","num-cds-protein","num-segs-target","num-introns","tot-cds-len","orig-seg-len","cds-div3","orig-seg-div3","orig-seg-sim",
    #             "cds-stops","seg-stops","stop-chunk-lens","stop-chunk-sims","final-seg-len","final-seg-sim","final-seg-div3"];
    # logfile.write("\t".join(headers) + "\n");
    # # The list of tracker dict keys to write to the file.
    # skippers = ["pid","sample","target","seg-start","seg-end","aa-seq","nt-seq","mm-seq","nt-aln","aa-aln","mm-aln"];
    #  # The list of tracker dict keys to NOT write to the file.
    

    # for pid in tracker_main:
    #     for target in tracker_main[pid]:
    #         for sample in tracker_main[pid][target]:
    #             outline = [pid, sample, target];
    #             # print(pid, sample, target);
    #             # print(tracker_main[pid][target][sample]);
    #             outline += [ str(tracker_main[pid][target][sample][value]) for value in headers if value not in skippers ];
    #             core.PWS("\t".join(outline), logfile);
    # # Write the number of CDS for each sample










'''
        for target in mus_pids[pid]['targets']:
            exonerate_file = os.path.join(indir, target + "-" + pid + "-exonerate.txt");
            seq_out_file = os.path.join(outdir, pid + ".fa");
            samples_found = {};
            target_tracker = {};

            cds_flag = False;
            seq_flag = False;
            for line in open(exonerate_file):
                print(line);
                if line == "# --- START OF GFF DUMP ---\n":
                    cds_flag = False;
                    seq_flag = False;

                if any(line.startswith(skip_str) for skip_str in ["Command line", "Hostname", "#", "--"]):
                    continue;
                # Skip lines that start with some key characters

                if not cds_flag:              
                    line = line.strip().split("\t");
                    if line[2] != "cds":
                        cds_flag = True;
                        cur_sample = line[0];

                        if cur_sample not in samples_found:
                            samples_found[cur_sample] = { 'starts' : [], 'ends' : [], 'lens' : [], 'seqs' : [] };

                if cds_flag:
                    split_line = line.strip().split("\t");
                    if split_line[0] == "INFO":
                        
                        if split_line[7] == "+":
                            start, end, length = int(split_line[8]), int(split_line[9]), int(split_line[10]);
                        elif split_line[7] == "-":
                            start, end, length = int(split_line[9]), int(split_line[8]), int(split_line[10]);

                        samples_found[cur_sample]['starts'].append(start);
                        samples_found[cur_sample]['ends'].append(end);
                        samples_found[cur_sample]['lens'].append(length);
                        cur_seq = "";

                    if line == "BEGIN SEQUENC":
                        seq_flag = True;
                        continue;

                    if seq_flag:
                        if line == "END SEQUENCE":
                            samples_found[cur_sample]['seqs'].append(cur_seq);
                            continue;
                        cur_seqs += line.strip();
                        
            for s in samples_found:
                print("SAMPLES FOUND")
                print(s, samples_found[s]);              
            return
'''







'''
    for gff_file in file_list:
    # Loop through all files
        if not gff_file.endswith(".gff"):
            skipped += 1;
            continue;
        # Skip any files that aren't GFF

        pid = gff_file[:gff_file.index("-")];
        mtid = mus_pid_to_tid[pid];
        if mtid not in mexp:
            skipped += 1;
            continue;
        rpid = mus_pid_to_rat_pid[pid];
        rtid = rat_pid_to_tid[rpid];
        tracker[pid] = {};
        # Get the protein ID from the file name and add it to the tracker

        gff_path = os.path.join(indir, gff_file);
        # Get the full path of the GFF file

        for line in open(gff_path):
        # Loop through all lines in current file
            if any(line.startswith(skip_str) for skip_str in ["Command line", "Hostname", "#", "--"]):
                continue;
            # Skip lines that start with some key characters

            line = line.strip().split("\t");
            try:
                if line[2] != "cds":
                    continue;
                sample, start, end, strand = line[0], line[3], line[4], line[6];
            except:
                print(" * ", gff_file, "\t", line);
                continue;
            # Parse the line and skip if not CDS

            if sample not in tracker[pid]:
                tracker[pid][sample] = { 'num-cds' : 0, 'starts' : [], 'ends' : [], 'seqs' : [], 'tot-len' : 0 };
            # Add the sample to the tracker for the current protein

            tracker[pid][sample]['starts'].append(start);
            tracker[pid][sample]['ends'].append(end);
            tracker[pid][sample]['num-cds'] += 1;
            tracker[pid][sample]['tot-len'] += (int(end) - int(start) + 1);
            # Add the CDS prediction info to the tracker for the current sample

            seqfile = os.path.join(seqdir, pid, sample + ".fa");
            # Use the sequence directory and the protein ID to get the original sequence file name for the current protein and sequence

            cmd = "samtools faidx " + seqfile + " " + sample + ":" + start + "-" + end
            if strand == "-":
                cmd += " -i";     
            cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
            assert cmd_result.stderr.decode() == "", "\nsamtools returned an error with the following command:\n" + cmd + "\n"
            # Retrieve the sequence of the current CDS with samtools

            cmd_out = cmd_result.stdout.decode()
            tracker[pid][sample]['seqs'].append("".join(list(filter(None, cmd_out.split("\n")[1:]))));
            # Parse the sequence from the samtools output and add it to the tracker

        outfilename = os.path.join(outdir, pid + ".fa");
        # Define the output file name for the current protein ID
        with open(outfilename, "w") as outfile:
            outfile.write(">mm10\n" + mseq[mtid] + "\n");
            add_rat = False;
            if rtid in rseq:
                add_rat = True;
                outfile.write(">rnor6\n" + rseq[rtid] + "\n");

            for sample in tracker[pid]:
            # Loop through every sample in the current protein ID

                sample_seq = "";
                # Initialize an empty sequence string

                sorted_indices = [ tracker[pid][sample]['starts'].index(start) for start in sorted(tracker[pid][sample]['starts'] ) ];
                # Sort the sequences for the current sample by sorting the starting coordinates and saving the sorted indices which will be used
                # to look up the sequences in the sorted order

                # print(tracker[pid][sample]);
                # print(sorted_indices);

                for ind in sorted_indices:
                    # print(tracker[pid][sample]['seqs'][ind]);
                    sample_seq += tracker[pid][sample]['seqs'][ind];
                # Loop through the sorted indices and add that sequence to the sequence string

                outfile.write(">" + sample + "\n" + sample_seq + "\n");
                # Write the sequence to the file

                tracker[pid][sample]['seqs'] = "";
                # Overwrite the sequences in the tracker with an empty string to free up memory

            tracker[pid]["mm10"] = {'starts' : [], 'ends' : [], 'seqs' : "", 'num-cds' : str(mexp[mtid]['num-cds']), 'tot-len' : str(mexp[mtid]['tot-len']) };
            if add_rat:
                tracker[pid]["rnor6"] = {'starts' : [], 'ends' : [], 'seqs' : "", 'num-cds' : str(rexp[rtid]['num-cds']), 'tot-len' : str(rexp[rtid]['tot-len']) };
            written += 1;

    return tracker, written, skipped;











                        # targ_strand_flag = True;  
                    # if target_strand == "-":
                    #     targ_strand_flag = False;
                    #     continue;         

                    # if cur_sample not in target_tracker[pid][target]:
                    #     target_tracker[pid][target][cur_sample] = { 'num-cds-protein' : 0, 'num-segs-target' : 0, 'num-introns' : [], 'seg-start' : [], 'seg-end' : [], 
                    #     'tot-cds-len' : 0,
                    #                                                  'orig-seg-len' : [], 'cds-div3' : "NA", 'orig-seg-div3' : "NA", 'orig-seg-sim' : [], 'cds-stops' : 0, 
                    #                                                  'seg-stops' : 0, 
                    #                                                  'nt-seq' : [], 'aa-seq' : [], 'mm-nt-seq' : [], 'mm-aa-seq' : [],
                    #                                                  'nt-aln' : [], 'aa-aln' : [], 'mm-nt-aln' : [], 'mm-aa-aln' : [],  
                    #                                                   'stop-chunk-sims' : "NA", 'stop-chunk-lens' : "NA",
                    #                                                   'final-seg-len' : "NA", 'final-seg-sim' : "NA", "final-seg-div3" : "NA" };
                    #                                                   #'nt-seq' : [], 'aa-seq' : [], 'mm-seq' : [],
                    #     # The main tracker dict tracks a lot. For each target/sample pair it retrieves sequences, alignments, coordinates,
                    #     # similarities, etc. If more than one segement is predicted for a target/sample pair it keeps all of them in lists. 
                    # # Initialize the current sample for the current target if it hasn't been already.

                    # if cur_sample not in sample_tracker:
                    #     sample_tracker[cur_sample] = {'nt-seq' : "", 'aa-seq' : "", "targets" : []};
                    #     # The sample_tracker dict only keeps track of sequences and number of targets with a predicted sequence (later
                    #     # counted as CDS exons). The sequences are concatenated for each target.
                    # # Initiize the sample dict for the current sample if it hasn't been already.
'''
