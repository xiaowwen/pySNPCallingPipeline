import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import utilities
from   read_config import read_config
from   parse_vcf import parse_vcf

def obtainSNPMatrix(config_file):

    config_reader      = read_config(config_file)
    config_dict        = config_reader.get_config_contents()
    reference          = config_dict['REF_FILE']
    sequences          = config_dict['INPUT_DIR']
    alignment_dir      = config_dict['ALIGNMENT_DIR']
    snp_checker_out    = config_dict['SNP_CHECKER_OUTPUT']
    snp_align_fasta    = config_dict['CREATE_ALIGN_OUTPUT']
    final_snp_out      = config_dict['CREATE_ALIGN_LIST_OUTPUT']

    unique_prefix_seqs = utilities.isolate_folders(sequences)
    output_vcf = alignment_dir + final_snp_out.split('/')[-1].split('.')[0] + '.vcf'
    snps_vcf   = parse_vcf(output_vcf)

    FINAL_SNP_POSITIONS = snps_vcf.VCF_DataFrame.position.values.tolist() # all positions from VCF file

    # Load all isolates, get their VCF objects and put it in a list here so that file IO doesn't happen at every step

    all_isolate_dfs = []
    for i in range(len(unique_prefix_seqs)):
        bwa_vcf_names  = alignment_dir + '/' + unique_prefix_seqs[i] + '/'  + unique_prefix_seqs[i] +"_sorted_bwa.vcf"
        vcf_reader     = parse_vcf(bwa_vcf_names)
        isolate_df_    = vcf_reader.VCF_DataFrame
        all_isolate_dfs.append([bwa_vcf_names,isolate_df_])

    # Finally loop over all reviewed positions and create the N_isolates by N_FINAL_SNP_POSITIONS matrix representing the SNP matrix

    SNPs_matrix = np.zeros([ len(all_isolate_dfs), len(FINAL_SNP_POSITIONS) ])
    isolate_list = []

    for j in range(len(FINAL_SNP_POSITIONS)):
        for k in range( len(all_isolate_dfs) ):
            isolate_SNPs = all_isolate_dfs[k][1]
            if FINAL_SNP_POSITIONS[j] in isolate_SNPs.position.values.tolist():
                SNPs_matrix[k,j] = 1

    for ij in range( len(all_isolate_dfs) ):
        isolate_path = all_isolate_dfs[ij][0]
        isolate_list.append(isolate_path.split('/')[-2])
    
    # Save SNP matrix and corresponding Isolate list to file

    SNP_pos_matrix = np.array(FINAL_SNP_POSITIONS)
    np.savetxt('SNP_Matrix.txt', SNPs_matrix.astype(int), fmt='%i', delimiter=",")
    np.savetxt('isolate_list.txt', isolate_list, fmt='%s', delimiter=" ")


# Run the analysis
obtainSNPMatrix("conf.intra.txt")
