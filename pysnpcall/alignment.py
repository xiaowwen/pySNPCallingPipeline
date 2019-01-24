import os
import glob
from pathlib import Path
import utilities
from   read_config import read_config
from   parse_vcf import parse_vcf

class alignment(object):

    # Reference file as to be reference.fasta or reference.fa
    # Sequence files must be sequences (paired end) files (PATH_TO_SEQUNCES)/*P.fq.gz (coming from read trimming) or (PATH_TO_SEQUNCES)/*.fq.

    def __init__(self, config_file):

        # Set class vaiables
        self.config_file      = config_file
        self.config_reader    = read_config(self.config_file)
        self.config_dict      = self.config_reader.get_config_contents()
        self.bwa_bin          = self.config_dict['BWA_PATH'] + '/' + 'bwa'
        self.last_bin         = self.config_dict['LAST_PATH'] + '/' + 'lastal'
        self.lastdb           = self.config_dict['LAST_PATH'] + '/' + 'lastdb'
        self.last_pair_prob   = self.config_dict['LAST_PATH'] + '/' + 'last-pair-probs'
        self.last_maf         = self.config_dict['LAST_PATH'][:-3] + '/scripts/' + 'maf-convert'
        self.novo_bin         = self.config_dict['NOVOALIGN_PATH'] + '/' + 'novoalign'
        self.novo_index       = self.config_dict['NOVOALIGN_PATH'] + '/' + 'novoindex' 
        self.samtools         = self.config_dict['SAMTOOLS_PATH'] + '/' + 'samtools'
        self.bcftools         = self.config_dict['BCFTOOLS_PATH'] + '/' +  'bcftools'

    def do_alignments(self):

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

        # Files for LAST aligner
        last_ref_dir = ref_dir + "last_ref/"
        last_reference = last_ref_dir + ref_fasta_prefix

        # Files for BWA alignner
        bwa_ref_dir = ref_dir + "bwa_ref/"
        bwa_reference = bwa_ref_dir + ref_fasta_prefix

        # Files for novoalign

        novo_ref_dir = ref_dir + "novo_ref/"
        novo_reference = novo_ref_dir + ref_fasta_prefix
        novo_index = ref_fasta_fasta + ".ndx"
        novo_index_path = novo_ref_dir + novo_index

        # Make reference/index directories 
        # LAST aligner
        os.system("mkdir %s; mkdir %s; cd %s; cp %s .; %s %s %s" % (ref_dir, last_ref_dir, last_ref_dir, reference, self.lastdb, ref_fasta_prefix, reference))
        # Novoalign 
        os.system("mkdir %s; cd %s; cp %s .; %s %s %s" % (novo_ref_dir, novo_ref_dir, reference, self.novo_index, novo_index, ref_fasta_fasta))

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
         
            bwa_vcf             = unique_prefix_seq + "_sorted_bwa.vcf"
            lastalign_vcf       = unique_prefix_seq + "_sorted_lastalign.vcf"
            novoalign_vcf       = unique_prefix_seq + "_sorted_novoalign.vcf"

            unique_prefix_dir   = alignment_dir + '/' + unique_prefix_seq

            if isolate_seq1.is_file() and isolate_seq2.is_file():

                #  1. LAST alignments
                ####################################################################################################

                print("Doing LAST alignments.....")
                os.system("mkdir %s" % (unique_prefix_dir))
                os.system("cd %s; %s -Q1 -e120 -i1 %s %s > rOne.maf" % (unique_prefix_dir, self.last_bin, last_reference, isolate_seq1))
                os.system("cd %s; %s -Q1 -e120 -i1 %s %s > rTwo.maf" % (unique_prefix_dir, self.last_bin, last_reference, isolate_seq2))
                os.system("cd %s; %s rOne.maf rTwo.maf > full_output.maf" % (unique_prefix_dir, self.last_pair_prob))
                os.system("cd %s; %s -d sam full_output.maf > full_output.sam" % (unique_prefix_dir, self.last_maf))

                # Convert SAM file into binary BAM and then to VCF

                os.system("cd %s; %s view -S -b full_output.sam > full_output.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s sort full_output.bam -o full_output_sorted.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s mpileup -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (unique_prefix_dir, self.bcftools, reference))
                os.system("cd %s; %s call --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (unique_prefix_dir, self.bcftools, lastalign_vcf))
                os.system("cd %s; rm rOne.maf rTwo.maf full_output.maf full_output.sam full_output.bam full_output_sorted.bam full_output_sorted.bcf" % (unique_prefix_dir))

                #  2. BWA alignments 
                ####################################################################################################

                print("Doing BWA alignments.....")
                os.system("cd %s; cp %s .; %s index %s" % (unique_prefix_dir, reference, self.bwa_bin, ref_fasta_fasta))
                os.system("cd %s; %s aln %s %s > rOne.sai" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq1))
                os.system("cd %s; %s aln %s %s > rTwo.sai" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq2))
                os.system("cd %s; %s sampe %s rOne.sai rTwo.sai %s %s > full_output.sam" % (unique_prefix_dir, self.bwa_bin, ref_fasta_fasta, isolate_seq1, isolate_seq2))

                # Convert SAM file into binary BAM and then to VCF

                os.system("cd %s; %s view -S -b full_output.sam > full_output.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s sort full_output.bam -o full_output_sorted.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s mpileup -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (unique_prefix_dir, self.bcftools, reference))
                os.system("cd %s; %s call --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (unique_prefix_dir, self.bcftools, bwa_vcf))
                os.system("cd %s; rm rOne.sai rTwo.sai full_output.sam full_output.bam full_output_sorted.bam full_output_sorted.bcf" % (unique_prefix_dir))

                #  3. Novoalign alignments 
                ####################################################################################################

                print("Doing Novoalign alignments.....")
                os.system("cd %s; cp %s %s .; gunzip *.gz" % (unique_prefix_dir, isolate_seq1, isolate_seq2))
                os.system("cd %s; %s -d %s -f %s %s -o SAM > full_output.sam" % (unique_prefix_dir, self.novo_bin, novo_index_path, isolate_seq1f, isolate_seq2f))

                # Convert SAM file into binary BAM and then to VCF

                os.system("cd %s; %s view -S -b full_output.sam > full_output.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s sort full_output.bam -o full_output_sorted.bam" % (unique_prefix_dir, self.samtools))
                os.system("cd %s; %s mpileup -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (unique_prefix_dir, self.bcftools, reference))
                os.system("cd %s; %s call --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (unique_prefix_dir, self.bcftools, novoalign_vcf))
                os.system("cd %s; rm full_output.sam full_output.bam *.fq full_output_sorted.bam full_output_sorted.bcf %s" % (unique_prefix_dir, index_files_for_removal))

#align = alignment('conf.txt')
#align.do_alignments()
