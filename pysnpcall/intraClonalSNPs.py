import pandas as pd
import numpy as np
from   read_config import read_config
from   parse_vcf import parse_vcf
import glob
import utilities

def intraClonalSNPs(config_file):

    # Read configuration file
    config_reader = read_config(config_file)
    config_dict   = config_reader.get_config_contents()

    # Obtain neccesary file paths, handles and parameters
    reference              = config_dict['REF_FILE']
    sequences              = config_dict['INPUT_DIR']
    alignment_dir          = config_dict['ALIGNMENT_DIR']
    intra_snp_dir          = config_dict['INTRA_SNP_WORKING_DIR']
    intra_snp_input        = config_dict['INTRA_SNP_INPUT']
    intra_snp_output       = config_dict['INTRA_SNP_OUTPUT']

    isolates = utilities.isolate_folders(sequences)
   
    # Get HQ SNPs list
    HQ_SNP_list = [int(snp.strip()) for snp in open(intra_snp_input)]
    HQ_SNP_list.sort()

    # Load ALL raw SNP calls
    all_isolate_dfs = []
    for i in range(len(isolates)):
        bwa_vcf_names  = alignment_dir + '/' + isolates[i] + '/'  + isolates[i] +"_sorted_bwa.vcf"
        vcf_reader     = parse_vcf(bwa_vcf_names)
        isolate_df_    = vcf_reader.VCF_DataFrame
        all_isolate_dfs.append(isolate_df_)

    intraclonal_HQ_SNPs = []
    # Circle through all VCFs and make a count of each time the HQ SNP occurs

    for HQ_SNP in HQ_SNP_list:

        count = 0
        for each_isolate_df in all_isolate_dfs:
            positions_in_df = each_isolate_df.position.values.tolist()
            
            if str(HQ_SNP) in positions_in_df:
                count += 1
        if count != len(all_isolate_dfs) and (count != 0):
            intraclonal_HQ_SNPs.append(HQ_SNP)

    # Write intra clonal SNPs to file

    SNP_pos_matrix = np.array(intraclonal_HQ_SNPs)
    np.savetxt(intra_snp_output, SNP_pos_matrix.astype(int), fmt='%i', delimiter=",")

    return intraclonal_HQ_SNPs

#intraClonalSNPs('conf.intra.txt')

