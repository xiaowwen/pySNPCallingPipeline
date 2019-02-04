import os
import glob
from pathlib import Path
import utilities
from   read_config import read_config
from   parse_vcf import parse_vcf

class getINDELs(object):

    # Reference file as to be reference.fasta or reference.fa
    # Sequence files must be sequences (paired end) files (PATH_TO_SEQUNCES)/*P.fq.gz (coming from read trimming) or (PATH_TO_SEQUNCES)/*.fq.

    def __init__(self, config_file):

        # Set class vaiables
        self.config_file      = config_file
        self.config_reader    = read_config(self.config_file)
        self.config_dict      = self.config_reader.get_config_contents()
        self.bwa_bin          = self.config_dict['BWA_PATH'] + '/' + 'bwa'
        self.samtools         = self.config_dict['SAMTOOLS_PATH'] + '/' + 'samtools'
        self.bcftools         = self.config_dict['BCFTOOLS_PATH'] + '/' +  'bcftools'
        self.gatk             = '/home/d/dguttman/emmanuel/Software/gatk-4.0.12.0/gatk' #'/home/enaziga/bin/gatk-4.0.6.0/gatk'
        self.picard           = '/home/d/dguttman/diazcaba/Software/picard.jar' #'/home/enaziga/bin/picard.jar'

    def do_getINDELs(self):

        # Obtain neccesary file paths, handles and parameters
        reference              = self.config_dict['REF_FILE']
        sequences              = self.config_dict['INPUT_DIR']
        alignment_dir          = self.config_dict['ALIGNMENT_DIR']

        # Prepare reference files
        ref_dir  = alignment_dir + "reference/"
        ref_file = reference.split('/')[-1]

        if ref_file.endswith('.fasta'):
            ref_fasta_prefix = ref_file[:-6]
        if ref_file.endswith('.fa'):
            ref_fasta_prefix = ref_file[:-3]

        ref_fasta_fasta  = ref_file

        # Files for BWA alignner
        bwa_ref_dir = ref_dir + "bwa_ref/"
        bwa_reference = bwa_ref_dir + ref_fasta_prefix

        # Grab sequence files
        try:
            files = glob.glob(sequences+"/*P.fq.gz")
        except:
            files = glob.glob(sequences+"/*.fq")
        files.sort()

        # Define text handles
        suffix1 = "_1P.fq.gz"
        suffix2 = "_2P.fq.gz"
        sf_novo1 = "_1P.fq"
        sf_novo2 = "_2P.fq"
        index_files_for_removal = ref_fasta_prefix+'.fa*'

        # Get file prefixes (Get paired end sequence files)

        unique_prefix_seqs = utilities.isolate_folders(sequences)
        
        for unique_prefix_seq in unique_prefix_seqs:

            print('Running Alignments for Replicon: ', unique_prefix_seq)

            R1             = sequences + '/' + unique_prefix_seq + suffix1
            R2             = sequences + '/' + unique_prefix_seq + suffix2
            R3             = unique_prefix_seq + sf_novo1
            R4             = unique_prefix_seq + sf_novo2
            isolate_seq1   = Path(R1)
            isolate_seq2   = Path(R2)
            isolate_seq1f  = Path(R3)
            isolate_seq2f  = Path(R4)
         
            bwa_vcf            = unique_prefix_seq + "_gatk_bwa.vcf"
            indels_vcf         = unique_prefix_seq + "_gatk_bwa_indels.vcf"
            filtered_vcf       = unique_prefix_seq + "_gatk_bwa_filtered.vcf"
            
            unique_prefix_dir  = alignment_dir + '/' + unique_prefix_seq
            isolate_slurm_file = unique_prefix_dir + '/' + unique_prefix_seq + "_get_indels.sh"

            if isolate_seq1.is_file() and isolate_seq2.is_file():
                
                os.system("mkdir %s" % unique_prefix_dir)
                #  2. BWA alignments 
                ####################################################################################################

                print("Carrying out BWA alignments.")
                os.system("cd %s; cp %s ." % (unique_prefix_dir, reference))
                string1 = "%s index %s" % (self.bwa_bin, ref_fasta_fasta)
                string2 = "%s aln %s %s > rOne.sai" % (self.bwa_bin, ref_fasta_fasta, isolate_seq1)
                string3 = "%s aln %s %s > rTwo.sai" % (self.bwa_bin, ref_fasta_fasta, isolate_seq2)
                string4 = "%s sampe %s rOne.sai rTwo.sai %s %s > full_output.sam" % (self.bwa_bin, ref_fasta_fasta, isolate_seq1, isolate_seq2)

                # Convert SAM file into binary BAM and then to VCF

                string5 = "%s view -S -b full_output.sam > full_output.bam" % (self.samtools)
                string6 = "%s sort full_output.bam -o full_output_sorted.bam" % (self.samtools)
                string7 = "%s mpileup -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (self.bcftools, ref_fasta_fasta)

                # Mark duplicates using PicardTools MarkDuplicates and call indels using GATK HaplotypeCaller
                
                string8 = "java -jar %s AddOrReplaceReadGroups I=full_output_sorted.bam O=AoR_full_output_sorted.bam RGID=HHMW3BGXY.1 RGPU=HHMW3BGXY.1.634 RGSM=634 RGPL=ILLUMINA RGLB=634 CREATE_INDEX=True" % (self.picard) 
                string9 = "java -jar %s MarkDuplicates I=AoR_full_output_sorted.bam O=duplicates.bam M=marked_dup_metrics.txt" % (self.picard) 

                #Run HaplotypeCaller (pao1 fasta was reference; indexed with samtools faidx v1.5 and Picard CreateSequenceDictionary v2.17.3; GATK version 3.8):
                string10 = "java -jar %s CreateSequenceDictionary R=%s " % (self.picard, ref_fasta_fasta)
                string11 = "%s faidx %s" % (self.samtools, ref_fasta_fasta)
                string12 = "%s index duplicates.bam" % (self.samtools)
                string13 = "%s --java-options \"-Xmx8G\"  HaplotypeCaller -R %s -I duplicates.bam -O %s -ploidy 1 -stand-call-conf 20" %(self.gatk, ref_fasta_fasta, bwa_vcf)
                string14 = "%s --java-options \"-Xmx8G\" SelectVariants -R %s -V %s -select-type INDEL -O %s" % (self.gatk, ref_fasta_fasta, bwa_vcf, indels_vcf) 

                string15 = "%s filter -i 'QUAL>=30 & DP>=20' %s -o %s " % (self.bcftools, indels_vcf, filtered_vcf) 
                string16 = "rm -f *.bam *.bcf *.sai *.sam" 

                all_strings = [string1,string2,string3,string4,string5,string6,string7,string8,string9,string10,string11,string12,string13,string14,string15,string16]

                # Write a slurm files for this isolate

                slurm1 = '#!/bin/bash'
                slurm2 = '#SBATCH --nodes=1'
                slurm3 = '#SBATCH --cpus-per-task=40'
                slurm4 = '#SBATCH --time=04:00:00'
                slurm5 = '#SBATCH --job-name=indels'
                slurm6 = '#SBATCH --output=indels_outlog.txt'
                slurm7 = 'module load java'
                slurm8 = 'module load python/3.6.5'

                slurm_lines = [slurm1,slurm2,slurm3,slurm4,slurm5,slurm6,'\n',slurm7,slurm8]

                with open(isolate_slurm_file, 'w') as outfile:

                    for slurm_line in slurm_lines:
                        outfile.write(slurm_line + '\n')

                    outfile.write('\n')\

                    for each_string in all_strings:
                        outfile.write(each_string + '\n')
            print(unique_prefix_dir)
            os.system("cd %s; chmod u+x %s; sbatch %s" % (unique_prefix_dir,isolate_slurm_file,isolate_slurm_file))
                
indels = getINDELs('indels_conf.txt')
indels.do_getINDELs()
