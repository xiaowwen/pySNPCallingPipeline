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
        self.gatk             = '/home/enaziga/bin/gatk-4.0.6.0/gatk'
        self.picard           = '/home/enaziga/bin/picard.jar'

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
         
            bwa_vcf           = unique_prefix_seq + "_gatk_bwa.vcf"
            filtered_vcf      = unique_prefix_seq + "_gatk_bwa_filtered.vcf"
            unique_prefix_dir = alignment_dir + '/' + unique_prefix_seq

            if isolate_seq1.is_file() and isolate_seq2.is_file():

                #  2. BWA alignments 
                ####################################################################################################

                print("Carrying out BWA alignments.")
                os.system("cd %s; cp %s .; %s index %s" % (unique_prefix_dir, reference, self.bwa_bin, ref_fasta_fasta))
                os.system("cd %s; %s aln %s %s > rOne.sai" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq1))
                os.system("cd %s; %s aln %s %s > rTwo.sai" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq2))
                os.system("cd %s; %s sampe %s rOne.sai rTwo.sai %s %s > full_output.sam" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq1, isolate_seq2))

                # Convert SAM file into binary BAM and then to VCF

                os.system("cd %s; %s view -S -b full_output.sam > full_output.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s sort full_output.bam -o full_output_sorted.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s mpileup -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (unique_prefix_dir, self.bcftools, reference))

                # Mark duplicates using PicardTools MarkDuplicates and call indels using GATK HaplotypeCaller

                os.system("cd %s; java -jar %s MarkDuplicates I=full_output_sorted.bam O=duplicates.bam M=marked_dup_metrics.txt" % (unique_prefix_dir, self.picard) )
                os.system("cd %s; java -jar %s CreateSequenceDictionary R=%s " % (unique_prefix_dir, self.picard, ref_fasta_fasta))
                os.system("cd %s; %s faidx %s" % (unique_prefix_dir, self.samtools, ref_fasta_fasta))
                os.system("cd %s; %s --java-options \"-Xmx4g\" HaplotypeCaller -R %s -I duplicates.bam -O %s \"-stand-call-conf 20\" \"-ploidy 1\" " %(unique_prefix_dir, self.gatk, ref_fasta_fasta, bwa_vcf))
                os.system("cd %s; bcftools filter -i 'QUAL>=30 & DP>=20 & TYPE=\"indel\" ' %s %s " % (unique_prefix_dir,bwa_vcf, filtered_vcf) )
                #os.system("cd %s; rm rOne.sai rTwo.sai full_output.sam full_output.bam full_output_sorted.bam duplicates.bam full_output_sorted.bcf" % (unique_prefix_dir))


indels = getINDELs('conf.intra.txt')
indels.do_getINDELs()
