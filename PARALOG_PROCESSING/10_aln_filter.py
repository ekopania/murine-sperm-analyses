#!/usr/bin/python
############################################################
# OLD INFO
# For Penn genomes, 06.2020
# Takes a log file from a clipkit run on amino acid sequence
# and removes corresponding sites from codon alignment.
# Taken from Gregg's murinae seq scripts and modified by E. Kopania 08.2022
############################################################

import sys
sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/software/core/corelib/");
# Add Gregg's corelib dir to the path.
sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/GREGG_MURINE_SCRIPTS/");
# Add Gregg's murine lib to the path.

import os, core, mseq, argparse
from Bio.Seq import Seq

############################################################
# Functions

def ntToCodon(nt_seq):
# This function takes a string of sequence and parses it into a list of codons.
# String is required to be divisible by 3.

    assert len(nt_seq) % 3 == 0, "\nOUT OF FRAME NUCLEOTIDE SEQUENCE! " + str(len(nt_seq));
    # Check that sequence is in frame.

    codon_seq = [(nt_seq[i:i+3]) for i in range(0, len(nt_seq), 3)];
    # Get chunks of 3 characters into a list.

    return codon_seq;

#########################

def countNongapLength(codon_seqs, seq_filter):
# This function goes through every sequence in an alignment and calculates
# the average length of each sequence excluding gaps.

    len_sum, num_seqs, num_gappy_seqs, gappy_seqs = 0, 0, 0, [];
    for title in codon_seqs:
        num_seqs += 1;
        full_len = len(codon_seqs[title]);
        nogap_len = len( [ codon for codon in codon_seqs[title] if codon != "---" ] );
        len_sum += nogap_len;

        if 1 - (nogap_len / full_len) > seq_filter:
            num_gappy_seqs += 1;
            gappy_seqs.append(title);

    return len_sum / num_seqs, num_gappy_seqs, gappy_seqs;

#########################

def countUniqIdentSeqs(codon_seqs):
# This function goes through every sequence in an alignment and counts how 
# many sequences are unique or identical.

    uniq_seqs, ident_seqs, found = 0, 0, [];
    codon_seq_list = list(codon_seqs.values());
    for seq in codon_seq_list:
        if codon_seq_list.count(seq) == 1:
            uniq_seqs += 1;
        if codon_seq_list.count(seq) != 1 and seq not in found:
            ident_seqs += 1;
            found.append(seq);

    return uniq_seqs, ident_seqs;

#########################

def siteCount(codon_seqs, aln_len, informative_sites_flag):
# This function goes through every site in an alignment to check for invariant sites, gappy sites, 
# and stop codons.

    stop_codons = ["TAG", "TAA", "TAR", "TGA", "TRA"];
    skip_chars = ["---", "NNN"];

    invar_sites, variant_sites, gap_sites, stop_codon_count, high_gap_sites = 0, 0, 0, 0, 0;
    informative_sites = "NA";
    if informative_sites_flag:
        informative_sites = 0;

    codon_seq_list = list(codon_seqs.values());
    #print(codon_seqs);
    for i in range(aln_len):
        site = [];
        for j in range(len(codon_seq_list)):
            site.append(codon_seq_list[j][i]);

            if informative_sites_flag:
                allele_counts = { allele : site.count(allele) for allele in site if allele not in skip_chars };
                # Count the occurrence of each allele in the site

                if len(allele_counts) > 1:
                    variant_sites += 1;
                    # If there is more than one allele in the site, it is variable

                    multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
                    # Count the number of allele present in at least 2 species

                    if len(multi_allele_counts) >= 2:
                        informative_sites += 1;
                    # If 2 or more alleles are present in 2 or more species, this site is informative
                else:
                    invar_sites += 1;
            # Counting informative sites if specified

            else:
                if site.count(site[0]) == len(site):
                    invar_sites += 1;
                else:
                    variant_sites += 1;
            # Otherwise just variant sites

        num_gaps = site.count("---");
        if num_gaps > 1:
            gap_sites += 1;

            if 1 - (num_gaps / len(site)) > 0.2:
                high_gap_sites += 1;

        stop_codon_count += len( [ codon for codon in site if codon in stop_codons ] );

    stop_codon_samples = [];
    for sample in codon_seqs:
        for codon in codon_seqs[sample]:
            if codon in stop_codons:
                stop_codon_samples.append(sample);
                break;

    return invar_sites, variant_sites, informative_sites, gap_sites, stop_codon_count, high_gap_sites, stop_codon_samples;

#########################

def codonWindowFilter(codon_seqs, num_samples, cds_len, site_filter):
# This function takes a sliding window along a codon alignment and filters windows
# where 2 or more codons have 2 or more gaps.

    exclude_sites = [];
    i = 0;
    while i < cds_len-3:
        sample_exclude = 0;
        #print(i);

        for sample in codon_seqs:
            #print(f, sample);
            c1 = codon_seqs[sample][i];
            c2 = codon_seqs[sample][i+1];
            c3 = codon_seqs[sample][i+2]; 

            gappy_c = [ c for c in [c1, c2, c3] if c.count("-") >= 2 ];

            if len(gappy_c) >= 2:
                sample_exclude += 1;

        if sample_exclude / num_samples > site_filter:
            exclude_sites += [i, i+1, i+2];

        i += 1;

    #exclude_sites.reverse();
    #print(exclude_sites);
    #print(cds_len, len(exclude_sites));

    for sample in codon_seqs:
        codon_seqs[sample] = [ codon_seqs[sample][i] for i in range(cds_len) if i not in exclude_sites ];

    return codon_seqs, len(list(set(exclude_sites)));

#########################

def writeAln(seqs, out_file):
# This function writes sequences to a file in FASTA format.

    with open(out_file, "w") as outfile:
        for sample in seqs:
            outfile.write(sample + "\n");
            outfile.write(seqs[sample] + "\n");

############################################################
# Options

parser = argparse.ArgumentParser(description="Codon alignment check filter");
parser.add_argument("-i", dest="input", help="Directory of alignments.", default=False);
#parser.add_argument("-w", dest="wsize", help="Codon window size. Default: 3", type=int, default=3);
parser.add_argument("-f", dest="pres_filter", help="The previous presence filter used in exonerate. For filenames.", default=False);
parser.add_argument("-s", dest="seq_filter", help="The gappy sequence filter threshold. Sequences that are greater than this percent gappy will be removed. Default: 20", type=int, default=20);
parser.add_argument("-c", dest="site_filter", help="The codon threshold. Sites that have greater than this percent sequences fail the codon window filter will be removed. Default: 50", type=int, default=50);
parser.add_argument("-o", dest="output", help="Desired output directory for filtered CDS alignments and translated AA alignments.", default=False);
#parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
# parser.add_argument("-e", dest="expected", help="The expected number of species in each alignment file. Check for one-to-one alignments to only align sequences that retained all species after trimming.", default=False);
#parser.add_argument("--noncoding", dest="noncoding", help="Set this option to check non-coding data. Will not check for stop codons.", action="store_true", default=False);
#parser.add_argument("--protein", dest="protein", help="Set this option to check amino acid data.", action="store_true", default=False);

parser.add_argument("--informativesites", dest="informative_sites", help="Set this option to also count informative sites per alignmend (takes much longer)", action="store_true", default=False);
parser.add_argument("--count", dest="count_only", help="Set this option to just provide the log file with counts/stats. Will not filter or write new sequences", action="store_true", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory with aligned CDS sequences must be defined with -i.");
args.input = os.path.abspath(args.input);

if not args.pres_filter:
    sys.exit( " * Error 2: Please provide the previous presence filter (-p).");

if args.seq_filter < 0 or args.seq_filter > 100:
    sys.exit( " * Error 3: Sequence gap threshold (-g) must be between 0 and 100.");

if args.site_filter < 0 or args.site_filter > 100:
    sys.exit( " * Error 4: Codon window gap threshold (-c) must be between 0 and 100.");

if not args.count_only and not args.output:
    sys.exit( " * Error 5: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
if args.output[-1] == "/":
    args.output = args.output[:-1];
args.output += "-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter);
# Adjust the output directory name to add the filter threshols.
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 6: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");

nt_outdir = os.path.join(args.output, "nt");
aa_outdir = os.path.join(args.output, "aa");
if not os.path.isdir(args.output):
    os.system("mkdir " + args.output);
if not args.count_only:
    os.system("mkdir " + nt_outdir);
    os.system("mkdir " + aa_outdir);
# Error checking and output file parsing.

# if args.noncoding and args.protein:
#     sys.exit(" * Error 4: Please specify only one of --noncoding or --protein.");
# elif args.noncoding:
#     mode = "nt";
# elif args.protein:
#     mode = "aa";
# else:
#     mode = "codon";

# expected = False;
# if args.expected:
#     eflag = False;
#     try:
#         expected = int(args.expected);
#     except:
#         eflag = True;
#     if eflag or expected < 1:
#         sys.exit(" * Error 4: -e must be a positive integer.");
# IO option error checking

pad = 26
# Job vars

sample_file = os.path.join(args.output, "sample-stats-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
log_file = os.path.join(args.output, "aln-stats-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# Job files

rm_stop_file = os.path.join(args.output, "stop-codon-filtered-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
rm_gappy_file = os.path.join(args.output, "gappy-seqs-filtered-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
rm_protein_file = os.path.join(args.output, "gappy-proteins-f" + args.pres_filter + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# Filter files

##########################
# Reporting run-time info for records.

with open(log_file, "w") as logfile, open(sample_file, "w") as samplefile:
    core.runTime("# CDS alignment filter", logfile);
    core.PWS("# IO OPTIONS", logfile);
    core.PWS(core.spacedOut("# Input CDS directory:", pad) + args.input, logfile);
    core.PWS(core.spacedOut("# Sequence gappiness threshold:", pad) + str(args.seq_filter), logfile);
    core.PWS(core.spacedOut("# Codon site gappiness threshold:", pad) + str(args.site_filter), logfile);
    #core.PWS(core.spacedOut("# Input sequence type:", pad) + mode, logfile);
    #core.PWS(core.spacedOut("# Codon window size:", pad) + str(args.wsize), logfile);
    core.PWS(core.spacedOut("# Output directory:", pad) + args.output, logfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", logfile);
    if args.count_only:
        core.PWS(core.spacedOut("# --count set:", pad) + "Will not output sequences.", logfile);
    core.PWS(core.spacedOut("# Log file:", pad) + log_file, logfile);
    core.PWS("# ----------------", logfile);
    
    args.seq_filter = args.seq_filter / 100;
    args.site_filter = args.site_filter / 100;

########################## 
# Filtering CDS aligns

    core.PWS("# " + core.getDateTime() + " Beginning filter.", logfile);

    rm_stop_codons, rm_gappy, rm_protein_gappy = [], [], [];
    pre_samples, pre_proteins = 0, 0;
    post_samples, post_proteins = 0, 0;

    aln_stats = ["num seqs", "codon aln length", "avg nongap length", "uniq seqs", "ident seqs", "gappy seqs", "invariant sites", "variant sites", "informative sites", "stop codons", 
                "percent sites with gap", "gappy sites"];
    aln_headers = ["align"] + [ "pre " + s for s in aln_stats ];
    if not args.count_only:
        aln_headers += ["sites filtered"] + [ "post " + s for s in aln_stats ];
    core.PWS("\t".join(aln_headers), logfile);
    # The alignment global headers

    sample_headers = ["sample"] + ["num alns", "num gappy"];
    samplefile.write("\t".join(sample_headers) + "\n");
    # The sample global headers

    sample_stats = {};
    # The sample counts dict
    
    fa_files = [ f for f in os.listdir(args.input) if f.endswith(".fa") ];
    num_alns = len(fa_files);
    num_alns_str = str(num_alns);
    num_alns_len = len(num_alns_str);
    # Read align file names from input directory

    init_sites = "NA";
    if args.informative_sites:
        init_sites = 0;

    written, num_short, num_high_ident, num_gappy, aln_prem_stop, num_no_info, num_stoppy = 0.0,0.0,0.0,0.0,0.0,0.0,0.0;
    pre_total_col, pre_total_nongap_col, pre_invariant, pre_variant, pre_informative, = 0,0,0,0,init_sites;
    post_total_col, post_total_nongap_col, post_invariant, post_variant, post_informative, = 0,0,0,0,init_sites;
    # Some count variables for all aligns

    # spec_high = {};
    # The dictionary to keep track of species count variables

    first_aln, counter, skipped = True, 0, 0;
    # Loop tracking variables

    for f in fa_files:
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_alns_len:
                counter_str = "0" + counter_str;
            print ("> " + core.getDateTime() + " " + counter_str + " / " + num_alns_str);
        counter += 1;
        # Loop progress   

        pre_proteins += 1;

        #print(f);

        cur_out = { h : "NA" for h in aln_headers };
        cur_out["align"] = f;
        # Initialize the current output dictionary.

        cur_infile = os.path.join(args.input, f);
        if not args.count_only:
            cur_nt_outfile = os.path.join(nt_outdir, f.replace(".fa", ".filter.fa"));
            cur_aa_outfile = os.path.join(aa_outdir, f.replace(".fa", ".filter.fa"));
        # Get the current in and output files

        seqs_orig = mseq.fastaGetDict(cur_infile);
        seqs = { t : seqs_orig[t].upper() for t in seqs_orig };
        samples = list(seqs.keys());
        pre_samples += len(samples);
        # Read the sequences

        for sample in seqs:
            if sample not in sample_stats:
                sample_stats[sample] = { col : 0 for col in sample_headers if col != "sample" };
            sample_stats[sample]['num alns'] += 1;
        # Count the samples in the alignment in the main dict and initialize if it is the first time this sample is seen

        codons = {};
        for seq in seqs:    
            codons[seq] = ntToCodon(seqs[seq]);
        # Convert from NT strings to codon lists

        #print(codons);

        num_seqs = len(codons);
        if num_seqs < 3:
            rm_protein_gappy.append(f);
            continue;

        cur_out["pre num seqs"] = str(num_seqs);
        # Count the number of samples in the alignment

        cur_out["pre codon aln length"] = len(codons[samples[0]]);
        pre_total_col += cur_out["pre codon aln length"];
        # Count the overall length of the alignment

        cur_out["pre avg nongap length"], cur_out["pre gappy seqs"], gappy_seqs = countNongapLength(codons, args.seq_filter);
        pre_total_nongap_col += cur_out["pre avg nongap length"];

        # if f != "ENSMUSP00000074707-mafft-cds.fa":
        #     continue;

        # print(f);
        # print(gappy_seqs);
        # print(len(gappy_seqs));
        # print(len(codons));

        # # for seq in codons:
        # #     print(seq, codons[seq]);
        # #     print(len(codons[seq]));
        # #     print(codons[seq].count("---"));
        # sys.exit();

        cur_out["pre uniq seqs"], cur_out["pre ident seqs"] = countUniqIdentSeqs(seqs);

        cur_out["pre invariant sites"], cur_out["pre variant sites"], cur_out["pre informative sites"], cur_out["pre percent sites with gap"], cur_out["pre stop codons"], cur_out["pre gappy sites"], pre_stop_samples = siteCount(codons, cur_out["pre codon aln length"], args.informative_sites);
        pre_invariant += cur_out["pre invariant sites"];
        pre_variant += cur_out["pre variant sites"];
        if args.informative_sites:
            pre_informative += cur_out["pre informative sites"];

        for sample in gappy_seqs:
            rm_gappy.append(f + "\t" + sample);

        if not args.count_only:

            f_codons = { sample : codons[sample] for sample in codons if sample not in gappy_seqs };

            cur_out["post num seqs"] = len(f_codons);

            if cur_out["post num seqs"] < 3:
                for col in aln_headers:
                    if col in ["align", "num seqs"]:
                        continue

                    col = "post " + col;
                    cur_out[col] = "NA";

            else:
                f_codons, cur_out["sites filtered"] = codonWindowFilter(f_codons, num_seqs, cur_out["pre codon aln length"], args.site_filter);

                f_samples = list(f_codons.keys());

                cur_out["post codon aln length"] = len(f_codons[f_samples[0]]);
                post_total_col += cur_out["post codon aln length"];

                cur_out["post avg nongap length"], cur_out["post gappy seqs"], gappy_seqs = countNongapLength(f_codons, args.seq_filter);
                post_total_nongap_col += cur_out["post avg nongap length"];

                cur_out["post uniq seqs"], cur_out["post ident seqs"] = countUniqIdentSeqs(f_codons);

                cur_out["post invariant sites"], cur_out["post variant sites"], cur_out["post informative sites"], cur_out["post percent sites with gap"], cur_out["post stop codons"], cur_out["post gappy sites"], post_stop_samples = siteCount(f_codons, cur_out["post codon aln length"], args.informative_sites);
                post_invariant += cur_out["post invariant sites"];
                post_variant += cur_out["post variant sites"];
                if args.informative_sites:
                    post_informative += cur_out["post informative sites"];

                for sample in post_stop_samples:
                    rm_stop_codons.append(f + "\t" + sample);

                post_proteins += 1;
                post_samples += len(f_codons);

                f_seqs = { sample : "".join(f_codons[sample]) for sample in f_codons if sample not in post_stop_samples };

                writeAln(f_seqs, cur_nt_outfile);

                aas = { sample : str(Seq(f_seqs[sample]).translate()) for sample in f_seqs };
                writeAln(aas, cur_aa_outfile);

        outline = [ str(cur_out[col]) for col in aln_headers ];
        logfile.write("\t".join(outline) + "\n");

    #print(sample_stats);
    for sample in sample_stats:
        outline = [sample] + [ str(sample_stats[sample][col]) for col in sample_headers if col != "sample" ];
        samplefile.write("\t".join(outline) + "\n");
    core.PWS("# ----------------", logfile);

    core.PWS("# Pre-filter proteins: " + str(pre_proteins), logfile);
    core.PWS("# Pre-filter samples : " + str(pre_samples), logfile);

    with open(rm_stop_file, "w") as stopfile:
        for seq in rm_stop_codons:
            stopfile.write(seq + "\n");
    core.PWS("# Samples removed with premature stops: " + str(len(rm_stop_codons)), logfile);

    with open(rm_gappy_file, "w") as gappyfile:
        for seq in rm_gappy:
            gappyfile.write(seq + "\n");
    core.PWS("# Samples removed above gappy threshold: " + str(len(rm_gappy)), logfile);

    with open(rm_protein_file, "w") as proteinfile:
        for protein in rm_protein_gappy:
            proteinfile.write(protein + "\n");
    core.PWS("# Proteins removed because of too few sequences post-filter: " + str(len(rm_protein_gappy)), logfile);

    core.PWS("# Post-filter proteins: " + str(post_proteins), logfile);
    core.PWS("# Post-filter samples : " + str(post_samples), logfile);
'''

        # if expected and num_seqs != expected:
        #     print("# Expected number of species not found... skipping: " + cur_infile + "\n");
        #     skipped += 1;
        #     continue;

        # if first_aln:
        # # Initialize some things on the first alignment
        #     spec_headers_order = [];
        #     for title in seqs:
        #         short_title = title.split(" ")[0];
        #         headers.append(short_title + " gaps/Ns");
        #         headers.append(short_title + " percent gaps/Ns");
        #         headers.append(short_title + " high percent gaps/Ns");
        #         if mode == "codon":
        #             headers.append(short_title + " prem stop");
        #         spec_headers_order.append(short_title);
        #         # For output, the species headers should retain an order set here

        #         spec_high[short_title] = { 'high-gaps' : 0, 'prem-stop' : 0 };
        #         # Initialize the overall counts for each species 

        #     logfile.write("\t".join(headers) + "\n");
        #     first_aln = False;
            # Write the headers and set the first flag to false

        num_high_gap, num_prem_stop = 0.0, 0.0;
        codon_seqs, spec_out = {}, {};
        seq_prem_stop = 0;
        first_seq, prem_stop_flag, short_seq, high_ident = True, False, "FALSE", "FALSE";
        # Variables for the current alignment

        for title in seqs:
            short_title = title.split(" ")[0];
            spec_out[short_title] = { 'gaps' : 0.0, 'perc-gaps' : 0.0, 'high-gaps' : "FALSE", 'prem-stop' : "FALSE" };
            # Initialize the species variables for this alignment

            if first_seq:
                seq_len = float(len(seqs[title]));
                cur_out["Aln length"] = str(seq_len);

                cur_out["Short seq"] = "FALSE";
                if seq_len < 100:
                    num_short += 1;
                    cur_out["Short seq"] = "TRUE";
                # Get the length of the alignment

                if mode == "codon":
                    codon_len = float(len(seqs[title]) / 3);
                    cur_out["Codon length"] = str(codon_len);
                #last_codon_ind = codon_len - 2;
                # Get the total number of codons in the alignment
                first_seq = False;
            # Get the length of the alignment from the first seq only

            seqs[title] = seqs[title].replace("!", "N");
            # Replace MACSE's frameshift ! char with N

            num_gaps = float(seqs[title].count("-") + seqs[title].count("N"));
            perc_gaps = round((num_gaps / seq_len) * 100, 2);
            if perc_gaps > 20:
                cur_out[">20% gap seqs"] += 1.0;
            # Count the number of gappy/uncalled sites for this sequence

            if mode == "codon":
                stop, seqs[title] = coreseq.premStopCheck(seqs[title], allowlastcodon=True, rmlast=True);
                if stop:
                    num_prem_stop += 1.0;
                    prem_stop_flag = True;
            # Check for premature stop codons and remove last codon if it is a stop
        # End sequence loop

        cur_out["Percent >20% gap seqs"] = round((cur_out[">20% gap seqs"] / num_seqs) * 100, 2);
        cur_out["Percent >20% gap seqs high"] = "FALSE";
        if cur_out["Percent >20% gap seqs"] > 50:
            cur_out["Percent >20% gap seqs high"] = "TRUE";
            num_gappy += 1;            

        all_seqs = list(seqs.values());
        ident_seqs = [];
        uniq_seqs = [];
        for t in seqs:
            seq_count = all_seqs.count(seqs[t]);
            if seq_count > 1 and seqs[t] not in ident_seqs:
                for i in range(seq_count):
                    ident_seqs.append(seqs[t]);
            elif seq_count == 1:
                uniq_seqs.append(seqs[t]);

        cur_out["Num uniq seqs"] = len(uniq_seqs) + len((set(ident_seqs)));
        cur_out["Num ident seqs"] = len(ident_seqs)

        cur_out["Percent ident seqs"] = round((cur_out["Num ident seqs"] / num_seqs) * 100, 2);
        cur_out["Percent ident seqs high"] = "FALSE";
        if cur_out["Percent ident seqs"] > 50:
            cur_out["Percent ident seqs high"] = "TRUE";
            num_high_ident += 1;
        # Check for identical sequences. This has to be done after the loop above because it changes sequences.

        if mode == "codon":
            if prem_stop_flag:
                aln_prem_stop += 1;
            # If any of the sequences in this alignment had a premature stop, add it to the count here

            cur_out["Num seq premature stop codons"] = num_prem_stop;
            cur_out["Percent seq premature stop codons"] = round((num_prem_stop / num_seqs) * 100, 2);
            cur_out["Percent seq premature stop codons high"] = "FALSE";
            if cur_out["Percent seq premature stop codons"] > 20:
                num_stoppy += 1;
                cur_out["Percent seq premature stop codons high"] = "TRUE";
            # Check if a high percentage of sequences contain premature stop codons
            
        cur_out["No info sites"] = 0.0;
        for i in range(int(seq_len)):
            num_gaps = 0.0;
            for title in seqs:
                if seqs[title][i] in ["-", "N"]:
                    num_gaps += 1.0;
            if num_gaps == num_seqs:
                num_no_info += 1;
                cur_out["No info sites"] += 1.0;
        # Count the number of columns that are all gaps or Ns

        cur_out["Percent no info sites"] = round((cur_out["No info sites"] / seq_len) * 100, 2);
        cur_out["Percent no info sites high"] = "FALSE";
        if cur_out["Percent no info sites"] > 20:
            num_gappy += 1;
            cur_out["Percent no info sites high"] = "TRUE";
        # Check if the number of columns that are all gaps or Ns is high

        outline = [ str(cur_out[h]) for h in headers ];
        logfile.write("\t".join(outline) + "\n");
        #core.PWS(("\t".join(outline)), logfile);
        # Write the log output line with counts for the current alignment

        if not args.count_only and not prem_stop_flag and cur_out["Short seq"] == "FALSE" and cur_out["Percent ident seqs high"] == "FALSE" and cur_out["Percent no info sites high"] == "FALSE":
            with open(cur_outfile, "w") as outfile:
                for title in seqs_orig:
                    outfile.write(title + "\n");
                    outfile.write(seqs_orig[title] + "\n");
                written += 1;
        # Write the edited sequence to the output file if there are no premature stop codons

    core.PWS("# ----------------", logfile);
    core.PWS(core.spacedOut("# Total aligns", 55) + str(num_alns), logfile);
    core.PWS(core.spacedOut("# Files skipped: ", pad) + str(skipped), logfile);
    core.PWS(core.spacedOut("# Aligns written", 55) + str(written), logfile);
    core.PWS(core.spacedOut("# Aligns shorter than 100bp", 55) + str(num_short), logfile);
    core.PWS(core.spacedOut("# Aligns with >50% identical sequences", 55) + str(num_high_ident), logfile);
    core.PWS(core.spacedOut("# Aligns with >20% of gappy/Ny sites", 55) + str(num_gappy), logfile);
    core.PWS(core.spacedOut("# Aligns with at least one premature stop:", 55) + str(aln_prem_stop), logfile); 
    if not args.noncoding:
        core.PWS(core.spacedOut("# Aligns with >20% of seqs with premature stops", 55) + str(num_stoppy), logfile);           
    core.PWS("# ----------------", logfile);
    # Write overall summary data

    # if args.noncoding:
    #     spec_headers = "Spec\tNumber >20% gappy"
    # else:
    #     spec_headers = "Spec\tNumber >20% gappy\tNumber premature stop"
    
    # core.PWS(spec_headers, logfile);
    # for title in spec_high:
    #     if mode == "codon":
    #         outline = [title, str(spec_high[title]['high-gaps']), str(spec_high[title]['prem-stop'])];
    #     else:
    #         outline = [title, str(spec_high[title]['high-gaps'])];
    #     core.PWS("\t".join(outline), logfile);
    # core.PWS("# ----------------", logfile);
    # Write species summary data

'''
