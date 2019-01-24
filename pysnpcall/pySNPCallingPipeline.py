import  argparse
import  numpy as np
from    read_config import read_config
from    parse_vcf import parse_vcf
import  utilities
from    alignment import alignment
from    scinet_alignment import scinet_alignment
import  getHQSNPs
import  checkSNPs
import  filterSNPs
import  write_vcf
import  intraClonalSNPs
from    pathlib import Path

parser = argparse.ArgumentParser(description='Run pySNPCallingPipeline')

# Require a confoguration file with all options for the pySNPCallingPipeline run
parser.add_argument('-c', '--conf', dest='CONF_FILE', required=True, type=str, help='A configuration file is required to run pySNPCallingPipeline.')

# Run alignment only
parser.add_argument('--aln_only', action='store_true', help='Alignment only, default is FALSE.')

# Run all other analysis given a set of alignments
parser.add_argument('--no_aln', action='store_true', help='Run all analysis using pre-run alignments, default is FALSE.')

# Run locally or on a supercomputing cluster
parser.add_argument('--local', action='store_true',  help='Is this a LOCAL run or should SLURM files be created? Default is True.')

# If running a supercomputing cluster, should only slurm files be created or should the jobs be submitted  as well
parser.add_argument('--submit', action='store_true', default=False, help='If running a supercomputing cluster, should only slurm files be created or should the jobs be submitted  as well.')

# Run a subset of analysis
parser.add_argument('--subset', dest='SUBSET', default=None,  help='Provide a comma seperated list of a subset of "alignment, getHQSNPs, intraClonalSNPs, checkSNPs, filterSNPs"')

# Default: run all analysis locally
parser.add_argument('--def_run', action='store_true',  help='Default: run entire calculation locally.')

args   = parser.parse_args()

# Get variables 
CONF_FILE = args.CONF_FILE
ALN_ONLY  = args.aln_only
NO_ALN    = args.no_aln
RUN_ENV   = args.local
SUBSET    = args.SUBSET
DEF_RUN   = args.def_run
SUBMIT    = args.submit

print(RUN_ENV)

print("Selected run options:  ")
print("====================================================================================================================")
print("====================================================================================================================")
print("Configuration file:                         ", CONF_FILE)
print("Alignment only (--aln_only):                ", ALN_ONLY)
print("No alignment (--no_aln):                    ", NO_ALN)
print("Run locally (--local):                      ", RUN_ENV)
print("Run a subset of the analysis (--subset):    ", SUBSET)
print("Submit slurm jobs (--submit):               ", SUBMIT)
print("Default run (--def_run):                    ", DEF_RUN)
print("====================================================================================================================")
print("====================================================================================================================")

# Function to check if alignment files are present

def check_alignments(isolates):

    # 1. Check that all alignment files are in the appropriate folders
    for isolate in isolates:
        # Handles for VCF files for the three alignment methods
        bwa_vcf             =  Path(alignment_dir + '/' + isolate + '/'  + isolate + "_sorted_bwa.vcf")
        lastalign_vcf       =  Path(alignment_dir + '/' + isolate + '/'  + isolate +"_sorted_lastalign.vcf")
        novoalign_vcf       =  Path(alignment_dir + '/' + isolate + '/'  + isolate +"_sorted_novoalign.vcf")

        if bwa_vcf.is_file() and lastalign_vcf.is_file() and novoalign_vcf.is_file():
            checkAlignmentFiles = True
        else:
            checkAlignmentFiles = False
            print('You are likely missing VCF files in some folders including that for replicon %s' % isolate)
            break
    return checkAlignmentFiles

# Turn off default if some partial analysis is requested
if NO_ALN or ALN_ONLY or SUBSET is not None:
    DEF_RUN = False
else:
    DEF_RUN = True
    print('Default: Running alignment on current workstation and carrying out SNP analysis:')

# ALIGNMENT only run
if ALN_ONLY:
    
    print('This is an alignment only run of pySNPCAllingPipeline')
    if RUN_ENV:
        align = alignment(CONF_FILE)
        align.do_alignments()
    else:
        # Prepare SLURM files
        align = scinet_alignment(CONF_FILE)
        align.do_scinet_alignments(submit_JOB=SUBMIT)

# DEFAULT RUN: Alignmnet --> getHQSNPs --> checkSNPs --> FilterSNPs --> write VCF
if DEF_RUN:
    
    # 1. Do alignments
    if RUN_ENV:
        align = alignment(CONF_FILE)
        align.do_alignments()
    else:
        align = scinet_alignment(CONF_FILE)
        align.do_scinet_alignments(submit_JOB=SUBMIT)
 
    # 2. Get HQ SNPs using alignments
    getHQSNPs.getHQSNPs(CONF_FILE)

    # 3. Get intraclonal SNPs

    intraClonalSNPs.intraClonalSNPs(CONF_FILE)

    # 4. Check HQ SNPs
    checkSNPs.checkSNPs(CONF_FILE)

    # 5. Filter SNPs
    filterSNPs.filterSNPs(CONF_FILE)
    print('SNP calling is complete')

    # 6. Write out a VCF file of the filtered SNPs
    write_vcf.write_vcf(CONF_FILE)

    print('VCF file written to the alignment directory')

# USING PRE-ALIGNED data: getHQSNPs --> checkSNPs --> FilterSNPs --> write VCF
if NO_ALN:

    print('Carrying out SNP analysis on previously aligned sequences. Checking VCF files....')

    # Gather some information from the configuration file
    config_reader = read_config(CONF_FILE)
    config_dict   = config_reader.get_config_contents()
    sequences     = config_dict['INPUT_DIR']
    alignment_dir = config_dict['ALIGNMENT_DIR']
    isolates      = utilities.isolate_folders(sequences)

    alignments_files_present = check_alignments(isolates)

    if alignments_files_present:

        # 2. Get HQ SNPs using alignments
        getHQSNPs.getHQSNPs(CONF_FILE)

        # 3. Get intraclonal SNPs
        intraClonalSNPs.intraClonalSNPs(CONF_FILE)

        # 4. Check HQ SNPs
        checkSNPs.checkSNPs(CONF_FILE)

        # 5. Filter SNPs
        filterSNPs.filterSNPs(CONF_FILE)
        print('SNP calling is complete')

        # 6. Write out a VCF file of the filtered SNPs
        write_vcf.write_vcf(CONF_FILE)

        print('VCF file written to the alignment directory')

# RUNNING A SUBSET of the pipeline
if SUBSET is not None:

    analysis_list = SUBSET.split(',')
    analysis_list = [word.lower() for word in analysis_list]

    # Determine analysis to be done
    analysis_arr  = np.zeros(5,dtype=int)
    if 'alignment' in analysis_list:
        analysis_arr[0] = 1
    if 'gethqsnps' in analysis_list:
        analysis_arr[1] = 1
    if 'intraclonalsnps' in analysis_list:
        analysis_arr[2] = 1
    if 'checksnps' in analysis_list:
        analysis_arr[3] = 1
    if 'filtersnps' in analysis_list:
        analysis_arr[4] = 1

    full_analysis_list = ['alignment','gethqsnps','intraclonalsnps', 'checksnps','filtersnps']

    executable_list = []
    for ijk in range(len(analysis_arr)):
        if analysis_arr[ijk] == 1:
            executable_list.append(full_analysis_list[ijk])

    print(analysis_list)
    print(executable_list)

    # If alignment is not part of the subset, check if files exist

    if 'alignment' in analysis_list:

        for analysis in executable_list:

            if analysis == 'alignment':

                # 1. Do alignments
                if RUN_ENV:
                    align = alignment(CONF_FILE)
                    align.do_alignments()
                else:
                    align = scinet_alignment(CONF_FILE)
                    align.do_scinet_alignments(submit_JOB=SUBMIT)

            if analysis == 'gethqsnps':
    
                # 2. Get HQ SNPs using alignments
                getHQSNPs.getHQSNPs(CONF_FILE)

            if analysis == 'intraclonalsnps':

                # 3. Get intraclonal SNPs
                intraClonalSNPs.intraClonalSNPs(CONF_FILE)

            if analysis == 'checksnps':

                # 4. Check HQ SNPs
                checkSNPs.checkSNPs(CONF_FILE)

            if analysis == 'filtersnps':
    
                # 5. Filter SNPs
                filterSNPs.filterSNPs(CONF_FILE)
                print('SNP calling is complete')

        # 6. Write out a VCF file of the filtered SNPs
        write_vcf.write_vcf(CONF_FILE)

        print('VCF file written to the alignment directory')


    else:
        # Gather some information from the configuration file
        config_reader = read_config(CONF_FILE)
        config_dict   = config_reader.get_config_contents()
        sequences     = config_dict['INPUT_DIR']
        alignment_dir = config_dict['ALIGNMENT_DIR']
        isolates      = utilities.isolate_folders(sequences)
        alignments_files_present = check_alignments(isolates)

        if not alignments_files_present:
        
            print('You are likely missing VCF files in some folders. Provide the neccessary files or request an alignment step.')

        else:

            for analysis in executable_list:

                if analysis == 'gethqsnps':

                    # 2. Get HQ SNPs using alignments
                    getHQSNPs.getHQSNPs(CONF_FILE)

                if analysis == 'intraclonalsnps':

                    # 3. Get intraclonal SNPs
                    intraClonalSNPs.intraClonalSNPs(CONF_FILE)

                if analysis == 'checksnps':

                    # 4. Check HQ SNPs
                    checkSNPs.checkSNPs(CONF_FILE)

                if analysis == 'filtersnps':

                    # 5. Filter SNPs
                    filterSNPs.filterSNPs(CONF_FILE)
                    print('SNP calling is complete')

            # 6. Write out a VCF file of the filtered SNPs
            write_vcf.write_vcf(CONF_FILE)

            print('VCF file written to the alignment directory')

