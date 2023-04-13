#!/usr/bin/python3
###########################################################
# For rodent exomes, 12.2020
# Runs coding2genome model of exonerate on each concatenated
# exome target to the appropriate mouse transcript.
###########################################################

import sys
sys.path.append("/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/GREGG_MURINE_SCRIPTS");
# Add the repo's lib dir to the path.

import os, mcore, argparse, re
from collections import defaultdict

###########################################################

exonerate_cmd_template = '''exonerate --model protein2genome --showalignment FALSE --showvulgar FALSE --showtargetgff -q {query} -t {target} --ryo "INFO\\t%ti\\t%tS\\t%tcb\\t%tce\\t%tcl\\t%s\\t%r$pi\\t%ps\\t%g\\t%qi\\t%qS\\t%qab\\t%qae\\t%qal\\t\\nBEGIN SEQUENCE\\n%tcsEND SEQUENCE\\nBEGIN MOUSE\\n%qasEND MOUSE\\nBEGIN ALN\\n{aln_str}END ALN\\n" &> {output}'''
# AA query 

#exonerate_cmd_template = '''exonerate --model coding2coding --showalignment FALSE --showvulgar FALSE --showtargetgff -q {query} -t {target} --ryo "INFO\\t%qi\\t%qS\\t%qcb\\t%qce\\t%qcl\\t%s\\t%r$pi\\t%ps\\t%g\\t%ti\\t%tS\\t%tab\\t%tae\\t%tal\\t\\nBEGIN SEQUENCE\\n%qcsEND SEQUENCE\\nBEGIN MOUSE\\n%tcsEND MOUSE\\nBEGIN ALN\\n{aln_str}END ALN\\n" &> {output}'''
# c2c

###########################################################
# Options

parser = argparse.ArgumentParser(description="Exonerate command generator");
parser.add_argument("-i", dest="input", help="Directory of input directories of FASTA files.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to exonerate. Default: exonerate", default="exonerate");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-nodes", dest="nodes", help="SLURM --nodes option.", default="1");
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-cpus", dest="cpus", help="SLURM --cpus-per-task option.", type=int, default=1);
parser.add_argument("-mem", dest="mem", help="SLURM --mem option.", type=int, default=0);
# SLURM options

args = parser.parse_args();
print(args.input);
if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory must be defined with -i.");
args.input = os.path.abspath(args.input);

if not args.name:
    name = mcore.getRandStr();
else:
    name = args.name;

if not args.output:
    sys.exit( " * Error 2: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
# if args.outname:
#     name = os.path.basename(args.output);
# else:
#     args.output = args.output + "-" + name + "/";
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if not args.part:
    sys.exit( " * Error 4: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 5: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 6: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -mem must be a positive integer.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

output_file = os.path.join(cwd, "jobs", name + ".sh");
submit_file = os.path.join(cwd, "submit", name + ".sh");
logdir = os.path.join(args.output, "logs");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    mcore.runTime("#!/bin/bash\n# Exonerate command generator", outfile);
    mcore.PWS("# IO OPTIONS", outfile);
    mcore.PWS(mcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    if args.outname:
        mcore.PWS(mcore.spacedOut("# --outname:", pad) + "Using end of output directory path as job name.", outfile);
    if not args.name:
        mcore.PWS("# -n not specified --> Generating random string for job name", outfile);
    mcore.PWS(mcore.spacedOut("# Job name:", pad) + name, outfile);
    mcore.PWS(mcore.spacedOut("# Output directory:", pad) + args.output, outfile);
    if args.overwrite:
        mcore.PWS(mcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        mcore.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    mcore.PWS(mcore.spacedOut("# Logfile directory:", pad) + logdir, outfile);
    if not os.path.isdir(logdir):
        mcore.PWS("# Creating logfile directory.", outfile);
        os.system("mkdir " + logdir);
    mcore.PWS(mcore.spacedOut("# Job file:", pad) + output_file, outfile);
    mcore.PWS("# ----------", outfile);
    mcore.PWS("# SLURM OPTIONS", outfile);
    mcore.PWS(mcore.spacedOut("# Submit file:", pad) + submit_file, outfile);
    mcore.PWS(mcore.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    mcore.PWS(mcore.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    mcore.PWS(mcore.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    mcore.PWS(mcore.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    mcore.PWS("# ----------", outfile);
    mcore.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.

    mouse_aa_dir = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/repro_genes-selected-cds-aa-trimmed/"
    #After trimming and framing exons

    #skip_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/skip_exons.txt";
    skip_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/scripts/mm10-exons-w-stops.txt";
    seqfiles = os.listdir(mouse_aa_dir);
    #targets_file = "../02-Annotation-data/Mus-selected-sequences_metadata_samp-counts_ratids-TESTSET.csv";
    #targets_file = "/mnt/beegfs/gt156213e/murinae-seq/02-Annotation/Mus-selected-sequences_metadata_samp-counts_ratids.csv";
    targets_file = "/uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/repro_exons_full.txt"
    mcore.PWS("# " + mcore.getDateTime() + " Reading IDs: " + targets_file);
    targets, targ_to_pid, pid_to_targ = {}, {}, defaultdict(list);
    first = True;
    for line in open(targets_file):
        if first:
            first = False;
            continue;
        
        #line = line.strip().split(",");
        line = line.strip().split("\t");
	#external_gene_name      ensembl_gene_id ensembl_transcript_id   ensembl_peptide_id      ensembl_exon_id chromosome_name
        targ, egn, tid, pid, eid = line[0].replace("\"", "")+"_"+line[4].replace("\"", ""), line[0].replace("\"", ""), line[2].replace("\"", ""), line[3].replace("\"", ""), line[4].replace("\"", "");
	#targ, coding, tid, pid, eid, num_targs = line[0].replace("\"", ""), line[1].replace("\"", ""), line[9].replace("\"", ""), line[13].replace("\"", ""), line[17].replace("\"", ""), line[22];
        #if pid == "NA" or num_targs == "NA" or coding == "FALSE":
        if pid == "NA":
            continue;

        targets[targ] = { 'egn' : egn, 'tid' : tid, 'pid' : pid, 'eid' : eid, 'mstrand' : "NA", "seqs" : {} };
    mcore.PWS("# Total targets: " + str(len(targets)));
    mcore.PWS("# ----------------");
    # Read the list of mouse target IDs.

    #skip_pids = [ line.strip().split(",")[0] for line in open(skip_file, "r") ];
    #skip_eids = [ line.strip().split(",")[1] for line in open(skip_file, "r") ];
    skip_ids = [ line.strip() for line in open(skip_file, "r") ];
    #print(skip_ids[1:5])

    num_exons, last_exon = 0, "";
    #aln_str = "{%Pqs\\t%Pqa\\t%Pqb\\t%Pqe\\t%Pts\\t%Pta\\t%Ptb\\t%Pte\\t%Pn\\t%Ps\\t%Pl\\n}";
    # AA

    aln_str = "{%Pqs\\t%Pqa\\t%Pqb\\t%Pqe\\t%Pts\\t%Pta\\t%Ptb\\t%Pte\\t%Pn\\t%Ps\\t%Pl\\n}";
    for target_id in targets:
    #for target_file in os.listdir(args.input):
        target_file_path = os.path.join(args.input, targets[target_id]['egn'] + ".fa");
        #target_id = target_file.replace(".fa", "");
        #target_file_path = os.path.join(args.input, target_file);
        query_id = targets[target_id]['eid'];
        query_pid_eid = ",".join((targets[target_id]['pid'], targets[target_id]['eid']));
        if query_pid_eid in skip_ids:
            print("skip")
            num_exons += 1;
            continue;

        r = re.compile(query_id);
        file_matches = list(filter(r.match, seqfiles));
        #print(file_matches);
        assert len(file_matches) == 1, "\nNO MATCHING FILES FOUND: " + query_id;

        query_file_path = os.path.join(mouse_aa_dir, file_matches[0]);

        #query_file_path = os.path.abspath(os.path.join(mouse_seq_dir, query_id + "-aa-trimmed.fa"));

        #print(num_proteins, query_id, query_pid, target_id);
        if query_pid_eid != last_exon:
            print(num_exons, query_pid_eid);
            last_exon = query_pid_eid;
            num_exons += 1;

        #base_input = os.path.splitext(query_dir)[0];
        #query_file = os.path.abspath(os.path.join(mouse_target_dir, query_dir + ".fa"));
        #cur_indir = os.path.join(args.input, query_dir);
        cur_outfile = os.path.join(args.output, target_id + "-exonerate.txt");
        #cur_logfile = os.path.join(logdir, base_input + "-mafft.log");

        #exonerate_cmd = "exonerate --model coding2genome " + protfilename + " " + infilename + " --bestn 5 --showtargetgff --showvulgar no --showalignment no --verbose 0 &> " + outfilename;
        #exonerate_cmd = "exonerate --model protein2genome --showalignment FALSE --showvulgar FALSE --showtargetgff -q " + query_file_path + " -t " + target_file_path 
        #exonerate_cmd += " --ryo '%qi\t%qS\t%qab\t%qae\t%qal\t%ti\t%tS\t%tcb\t%tce\t%tcl\t%s\t%r$pi\t%ps\t%g\n%tcs\n%qas\n'"
        #exonerate_cmd += " &> " + cur_outfile;
        exonerate_cmd = exonerate_cmd_template.format(query=query_file_path, target=target_file_path, aln_str=aln_str, output=cur_outfile);

        #exonerate --model coding2genome --ryo ">%ti (%tab - %tae)\n%tas\n>%qi (%qab - %qae)\n%qas" --showalignment FALSE --showvulgar TRUE --softmasktarget TRUE --score 1000 --bestn 5 --querytype dna --targettype dna --query QUERIES --target targeted_txpts.fa

        outfile.write(exonerate_cmd + "\n");

        # if num_proteins >= 900:
        #    break;

##########################
# Generating the submit script.

with open(submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.kopania@genetics.utah.edu
#SBATCH --account=clarkn-rw
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={tasks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

module load exonerate

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=args.tasks, cpus=args.cpus, mem=args.mem, output_file=output_file));

##########################







































# ##########################
# # Generating the commands in the job file.

#     for spec in genomes:
#         num_skipped = 0;
#         print(" --> " + spec)
#         for f in os.listdir(genomes[spec]['dir']):
#             if not f.endswith(".fa"):
#                 continue;

#             tid = f[:f.index("-")];
#             infilename = os.path.join(genomes[spec]['dir'], f);
#             outfilename = os.path.join(genomes[spec]['dir'], tid + "-exonerate.gff");
#             protfilename = os.path.join(mm10dir, tid + "-aa.fa");

#             if not os.path.isfile(protfilename):
#                 num_skipped += 1;
#                 continue;

#             #assert os.path.isfile(protfilename), "\nCouldn't find protein file: " + protfilename;

#             exonerate_cmd = "exonerate --model protein2genome " + protfilename + " " + infilename + " --bestn 5 --showtargetgff --showvulgar no --showalignment no --verbose 0 &> " + outfilename;
#             outfile.write(exonerate_cmd + "\n");
#         mcore.PWS("# " + spec + " : protein not found : " + str(num_skipped), outfile);

# #exonerate --model protein2genome ENSMUST00000000001-aa.fa ENSMUST00000000001-gdol.fa --bestn 1 --showtargetgff --showvulgar no --showalignment no --verbose 0 > test.gff

# ##########################
# # Generating the submit script.

# with open(submit_file, "w") as sfile:
#     submit = '''#!/bin/bash
# #SBATCH --job-name={name}
# #SBATCH --output={name}-%j.out
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=gregg.thomas@umontana.edu
# #SBATCH --partition={partition}
# #SBATCH --nodes=1
# #SBATCH --ntasks={tasks}
# #SBATCH --cpus-per-task={cpus}
# #SBATCH --mem={mem}

# parallel -j {tasks} < {output_file}'''

#     sfile.write(submit.format(name=args.name, partition=args.part, tasks=args.tasks, cpus=args.cpus, mem=args.mem, output_file=output_file));

# ##########################
