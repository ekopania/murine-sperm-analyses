#get gene fastas that only include species with RTM data
get_RTM_samples.sh

#Trim mm10 reference exons and get them in frame so that exonerate works properly
time -p python 06_frame_exons.py > logs/frame_exons_debug.log

#Run exonerate
#time -p python 07_exonerate_gen_2.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/GENE_FASTAS/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE_OUTPUT -n exonerate_repro_paralogs -part good_lab_reincarnation -tasks 100 > logs/exonerate_gen.log
#Run exonerate on RTM samples ONLY
time -p python 07_exonerate_gen_2.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/GENE_FASTAS_RTM/ -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE_OUTPUT_RTM -n exonerate_repro_paralogs_RTM -part clarkn-shared-rw -tasks 60 > logs/exonerate_gen.RTM.log
#Run exonerate on BLAST filtered gene fastas
#time -p python 07_exonerate_gen_2.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/GENE_FASTAS_FILTERED/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE_OUTPUT_BLASST_FILTERED -n exonerate_repro_paralogs -part good_lab_reincarnation_blast_filtered -tasks 100 > logs/exonerate_gen.blast_filtered.log

#Count number of sample hits for each exon and make a new file with this info; needed for exonerate2cds
#EDIT THIS SCRIPT FOR DIFFERENT INPUTS
sbatch 07b_get_exon_counts.sh

#Parse exonerate output; -f 50 because 50 is median from 07b_get_exon_counts.sh
#time -p python 08_exonerate_to_cds_2_trimmed.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE_OUTPUT -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE2CDS -f 50 > logs/exonerate-to-cds-f50.log
#RTM only; -f 16 because 16 is median
#time -p python 08_exonerate_to_cds_2_trimmed.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE_OUTPUT_RTM -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM -f 16 > logs/exonerate-to-cds-f16.RTM.log
#RTM only; parse by EXON rather than by PROTEIN
time -p python 08b_exonerate_to_cds_2_trimmed.byExon.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE_OUTPUT_RTM -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_byExon_RTM -f 16 > logs/exonerate-to-cds-f16.byExon.RTM.log
#BLAST filtered version
#time -p python 08_exonerate_to_cds_2_trimmed.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE_OUTPUT_BLAST_FILTERED -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE2CDS_BLAST_FILTERED -f 175 > logs/exonerate-to-cds-f175.blast_filtered.log
#Parse exonerate output WITHOUT mm10 ref seq in output fasta
#time -p python 08b_exonerate_to_cds_2_trimmed.NOmm10.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE_OUTPUT -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE2CDS_NOmm10 -f 0 > logs/exonerate-to-cds-f0.NOmm10.log

#merge paralogs
#EDIT THIS SCRIPT FOR DIFFERENT INPUTS
#Also, run twice, once for nt and once for aa
python merge_paralogs.py

#Align with mafft
# mkdir ../MAFFT
time -p python mafft_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/aa-merge-paralogs/ -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/MAFFT_RTM-f16/aa/ -n exonerate_mafft_f16_repro_paralogs_RTM -part clarkn-shared-rw -tasks 51 > logs/mafft_gen.f16.RTM.log
#WITHOUT mm10
#time -p python mafft_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE2CDS_NOmm10-f0/aa/ -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/MAFFT_NOmm10/aa/ -n exonerate_mafft_repro_paralogs.NOmm10 -part good_lab_reincarnation -tasks 50 > logs/mafft_gen.NOmm10.log

#Backtranslate from aa to nt
time -p python 09_backtranslate.py -aa /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/MAFFT_RTM-f16/aa/ -nt /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/EXONERATE2CDS_RTM-f16/nt-merge-paralogs/ -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/MAFFT_RTM-f16/nt/ > logs/backtranslate.f16.RTM.log
#WITHOUT mm10
#time -p python 10_backtranslate.py -aa /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/MAFFT_NOmm10/aa/  -nt /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/EXONERATE2CDS_NOmm10-f0/nt/  -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/MAFFT_NOmm10/nt/ > logs/backtranslate.NOmm10.log

#Alignment filtering
#Defaults: seqs more than 20% gappy removed; removes almost all seqs for these paralogs; only 8 proteins remaining
time -p python 10_aln_filter.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/MAFFT_RTM-f16/nt/ -f 16 -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM > logs/aln_filter.f16.RTM.log
#Seqs more than 50% gappy removed; 18 proteins remaining
time -p python 10_aln_filter.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/MAFFT/nt/ -f 0 -s 50 -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED > logs/aln_filter_s50.log
#Seqs more than 75% gappy removed; 30 proteins remaining
time -p python 10_aln_filter.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/MAFFT/nt/ -f 0 -s 75 -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED > logs/aln_filter_s75.log

#Run this to remove paralogs - if multiple seqs for a given sample and exon, only keep the one with the highest exonerate score
python remove_paralogs.py > remove_paralogs.log

#Generate tree based on default filter
time -p python /uufs/chpc.utah.edu/common/HIPAA/u6035720/software/core/generators/iqtree_gt_gen.py -i /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/FILTERED_RTM-f16-seq20-site50/nt -o /uufs/chpc.utah.edu/common/HIPAA/u6035720/MURINAE/PARALOG_PROCESSING/IQTREE_RTM-f16-seq20-site50 -b 1000 -n iqtree_repro_paralogs_RTM-f16-seq20-site50 -part clarkn-shared-rw -tasks 12 -p /uufs/chpc.utah.edu/common/HIPAA/u6035720/software/iqtree-2.2.0-Linux/bin/iqtree2 > logs/iqtree_RTM-f16-seq20-site50.log
#Generate tree based on 50% gappy filter
time -p python /home/ek112884/software/core/generators/iqtree_gt_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED-f0-seq50-site50/nt -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE-f0-seq50-site50 -b 1000 -n iqtree_repro_paralogs-f0-seq50-site50 -part good_lab_reincarnation -tasks 18 > logs/iqtree-f0-seq50-site50.log
#Generate tree based on 75% gappy filter
time -p python /home/ek112884/software/core/generators/iqtree_gt_gen.py -i /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/FILTERED-f0-seq75-site50/nt -o /mnt/beegfs/ek112884/murinae/PARALOG_PROCESSING/IQTREE-f0-seq75-site50 -b 1000 -n iqtree_repro_paralogs-f0-seq75-site50 -part good_lab_reincarnation -tasks 30 > logs/iqtree-f0-seq75-site50.log

#Get rid of protein ID appended to ends of sample names in final filtered fasta - necessary for BS test generator to work properly
#Run from inside FILTERED/nt/ directory
ls *fa | while read file; do sed -i 's/_|_.*//' ${file}; done
#ls *fa | while read file; do sed -i 's/_ENSMUSP.*//' ${file}; done #For mm10
#Run from inside IQTREE/loci/ directory
ls *filter/codeml.tre | while read file; do sed -i 's/_|_Mus_musculus_[[:alnum:]]\+-p[[:digit:]]\+:/:/g' ${file}; done

