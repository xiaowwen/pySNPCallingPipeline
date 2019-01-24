import os
import glob
from pathlib import Path
import utilities
from   read_config import read_config
from   parse_vcf import parse_vcf

class scinet_alignment(object):

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
        self.jobs_dir         = self.config_dict['JOBS_DIR'] + '/'

    def do_scinet_alignments(self, submit_JOB=False):

        # Obtain neccesary file paths, handles and parameters
        reference              = self.config_dict['REF_FILE']
        sequences              = self.config_dict['INPUT_DIR']
        alignment_dir          = self.config_dict['ALIGNMENT_DIR']
        jobfile                = self.jobs_dir + 'slurm_jobs.sh'

        # Prepare reference files
        ref_dir  = alignment_dir + "reference/"
        ref_file = reference.split('/')[-1]

        if ref_file.endswith('.fasta'):
            ref_fasta_prefix = ref_file[:-6]
        if ref_file.endswith('.fa'):
            ref_fasta_prefix = ref_file[:-3]

        ref_fasta_fasta  = ref_file
        ref_fasta_fasta_star = ref_file + ".*"

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
        os.system("mkdir %s; cd %s; cp %s .; %s %s %s > err" % (novo_ref_dir, novo_ref_dir, reference, self.novo_index, novo_index, ref_fasta_fasta))

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

            unique_prefix_dir   = alignment_dir + '/' + unique_prefix_seq

            bwa_slurm_file      = unique_prefix_dir + '/' + unique_prefix_seq + "_slurm_bwa.sh"
            last_slurm_file     = unique_prefix_dir + '/' + unique_prefix_seq + "_slurm_last.sh"
            novo_slurm_file     = unique_prefix_dir + '/' + unique_prefix_seq + "_slurm_novo.sh"
            combined_slurm      = unique_prefix_dir + '/' + "combined_slurm_file.sh"
            
            bwa_vcf             = unique_prefix_dir + '/' + unique_prefix_seq + "_sorted_bwa.vcf"
            lastalign_vcf       = unique_prefix_dir + '/' + unique_prefix_seq + "_sorted_lastalign.vcf"
            novoalign_vcf       = unique_prefix_dir + '/' + unique_prefix_seq + "_sorted_novoalign.vcf"

            if isolate_seq1.is_file() and isolate_seq2.is_file():

                #  1. LAST alignments
                ####################################################################################################

                os.system("mkdir %s" % (unique_prefix_dir))

                string1 ="%s -P 20 -Q1 -e120 -i1 %s %s > rOne.maf" % (self.last_bin, last_reference, isolate_seq1)
                string2 ="%s -P 20 -Q1 -e120 -i1 %s %s > rTwo.maf" % (self.last_bin, last_reference, isolate_seq2)
                string3 ="%s rOne.maf rTwo.maf > full_output.maf" % (self.last_pair_prob)
                string4 ="%s -d sam full_output.maf > full_output.sam" % (self.last_maf)

                # Convert SAM file into binary BAM and then to VCF

                string5 ="%s view --threads 20 -S -b full_output.sam > full_output.bam" % (self.samtools)
                string6 ="%s sort --threads 20 full_output.bam -o full_output_sorted.bam" % (self.samtools)
                string7 ="%s mpileup --threads 20 -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (self.bcftools, reference)
                string8 ="%s call --threads 20 --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (self.bcftools, lastalign_vcf)
                string9 ="rm rOne.maf rTwo.maf full_output.maf full_output.sam full_output.bam full_output_sorted.bam full_output_sorted.bcf"# % (unique_prefix_dir)

                last_strings = [string1,string2,string3,string4,string5,string6,string7,string8,string9]
                
                with open(last_slurm_file, 'w') as outfile1:
                    for last_string in last_strings:
                        outfile1.write(last_string + '\n')


                #  2. BWA alignments 
                ####################################################################################################
                
                #os.system("cp %s %s; %s index %s" % (reference, unique_prefix_dir, self.bwa_bin, reference)) #ref_fasta_fasta))
                string0 ="cp %s %s" % (reference, unique_prefix_dir)
                stringX ="%s index %s" % (self.bwa_bin, ref_fasta_fasta)
                string1 ="%s aln %s %s > rOne.sai" % (self.bwa_bin, ref_fasta_fasta, isolate_seq1)
                string2 ="%s aln %s %s > rTwo.sai" % (self.bwa_bin, ref_fasta_fasta, isolate_seq2)
                string3 ="%s sampe %s rOne.sai rTwo.sai %s %s > full_output.sam" % (self.bwa_bin, ref_fasta_fasta, isolate_seq1, isolate_seq2)

                # Convert SAM file into binary BAM and then to VCF

                string4 ="%s view --threads 20 -S -b full_output.sam > full_output.bam" % (self.samtools)
                string5 ="%s sort --threads 20 full_output.bam -o full_output_sorted.bam" % (self.samtools)
                string6 ="%s mpileup --threads 20 -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (self.bcftools, reference)
                string7 ="%s call --threads 20 --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (self.bcftools, bwa_vcf)
                string8 ="rm rOne.sai rTwo.sai full_output.sam full_output.bam full_output_sorted.bam full_output_sorted.bcf"# % (unique_prefix_dir)

                bwa_strings = [string0, stringX, string1, string2,string3,string4,string5,string6,string7,string8]
                with open(bwa_slurm_file, 'w') as outfile2:
                    for bwa_string in bwa_strings:
                        outfile2.write(bwa_string + '\n')

                #  3. Novoalign alignments 
                ####################################################################################################

                os.system("cp %s %s %s; cd %s; gunzip *.gz" % ( isolate_seq1, isolate_seq2, unique_prefix_dir, unique_prefix_dir))

                #os.system("mkdir %s; cd %s; cp %s .; %s %s %s" % (novo_ref_dir, novo_ref_dir, reference, self.novo_index, novo_index, ref_fasta_fasta))

                string1 ="%s -d %s -f %s %s -o SAM > full_output.sam" % (self.novo_bin, novo_index_path, isolate_seq1f, isolate_seq2f)

                # Convert SAM file into binary BAM and then to VCF

                string2 ="%s view --threads 20 -S -b full_output.sam > full_output.bam" % (self.samtools)
                string3 ="%s sort --threads 20 full_output.bam -o full_output_sorted.bam" % (self.samtools)
                string4 ="%s mpileup --threads 20 -Ob -o full_output_sorted.bcf -f %s full_output_sorted.bam" % (self.bcftools, reference)
                string5 ="%s call --threads 20 --ploidy 1 -vmO v -o %s full_output_sorted.bcf" % (self.bcftools, novoalign_vcf)
                string6 ="rm full_output.sam full_output.bam *.fq full_output_sorted.bam full_output_sorted.bcf"# % (unique_prefix_dir)

                novo_strings = [string1, string2,string3,string4,string5,string6]
                with open(novo_slurm_file, 'w') as outfile3:
                    for novo_string in novo_strings:
                        outfile3.write(novo_string + '\n')

            # Write a slurm files for this isolate

            slurm1 = '#!/bin/bash'
            slurm2 = '#SBATCH --nodes=1'
            slurm3 = '#SBATCH --cpus-per-task=40'
            slurm4 = '#SBATCH --time=04:00:00'
            slurm5 = '#SBATCH --job-name=alignments'
            slurm6 = '#SBATCH --output=alignments_outlog.txt'
            slurm7 = 'module load java'
            slurm8 = 'module load python/3.6.5'
            slurm_lines = [slurm1,slurm2,slurm3,slurm4,slurm5,slurm6,'\n',slurm7,slurm8]

            with open(combined_slurm, 'w') as outfile4:

                for slurm_line in slurm_lines:
                    outfile4.write(slurm_line + '\n')

                all_slurm_files = [bwa_slurm_file, novo_slurm_file, last_slurm_file]
                for slurm_file in all_slurm_files:
                    slurm_string = ' ' + slurm_file
                    outfile4.write(slurm_string + '\n')
                clean_up_ref = "rm " + ref_fasta_fasta + "  " + ref_fasta_fasta_star
                outfile4.write(clean_up_ref + '\n')

            os.system("cd %s; chmod u+x *.sh" % unique_prefix_dir)

            if submit_JOB:
                os.system("cd %s; sbatch %s" % (unique_prefix_dir, combined_slurm) )


#align = scinet_alignment('conf.intra.txt')
#align.do_scinet_alignments(submit_JOB=True)
