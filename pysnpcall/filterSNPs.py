import pandas as pd
import numpy as np
from   read_config import read_config
from   parse_vcf import parse_vcf
import glob
import utilities

def filterSNPs(config_file):

    # Read configuration file
    config_reader = read_config(config_file)
    config_dict   = config_reader.get_config_contents()

    # Obtain neccesary file paths, handles and parameters
    reference              = config_dict['REF_FILE']
    sequences              = config_dict['INPUT_DIR']
    alignment_dir          = config_dict['ALIGNMENT_DIR']
    snp_filter_input       = config_dict['SNP_FILTER_INPUT']
    snp_checker_out        = config_dict['SNP_CHECKER_OUTPUT']
    snp_filter_output      = config_dict['SNP_FILTER_OUTPUT']
    snp_filter_quality     = float(config_dict['SNP_FILTER_MIN_GOOD_QUALITY'])
    snp_filter_call_ratio  = float(config_dict['SNP_FILTER_MIN_GOOD_CALL_RATIO'])
    snp_filter_balance     = config_dict['SNP_FILTER_MIN_GOOD_BALANCE']

    snp_align_input        = config_dict['CREATE_ALIGN_INPUT']
    snp_align_fasta        = config_dict['CREATE_ALIGN_OUTPUT']
    final_snp_out          = config_dict['CREATE_ALIGN_LIST_OUTPUT']

    isolates               = utilities.isolate_folders(sequences)
    number_of_isolates     = float(len(isolates))
    checked_SNPs           = pd.read_csv(snp_checker_out)
    checked_SNPs_for_rev   = checked_SNPs.values.tolist()

    reviewed_SNP_list    = []
    reviewed_SNP_strings = []
    positions_only       = []

    print("Number of SNPs before filtering: ", len(checked_SNPs) )

    # Review every HQ SNP by loading up the slice of the data containing that SNP from the checkedSNP dataframe
    for SNPs_to_be_checked in checked_SNPs_for_rev:
        #print(SNPs_to_be_checked)
        
        isolate    = SNPs_to_be_checked[0]
        position   = SNPs_to_be_checked[1]
        alt_base   = SNPs_to_be_checked[2]
        ref_base   = SNPs_to_be_checked[3]

        quality    = float(SNPs_to_be_checked[4])
        readDepth  = float(SNPs_to_be_checked[5])
        refForward = float(SNPs_to_be_checked[6])
        refReverse = float(SNPs_to_be_checked[7])
        altForward = float(SNPs_to_be_checked[8])
        altReverse = float(SNPs_to_be_checked[9])
        refRatio   = float(SNPs_to_be_checked[10])

        refTotal   = refForward + refReverse
        altTotal   = altForward + altReverse
        refRatio   = refTotal/(altTotal+refTotal)

        refRequiredRatio = snp_filter_call_ratio
        altRequiredRatio = 1.0 - snp_filter_call_ratio

        isolate_name     = isolate.split('/')[-2]
    
        # Reference and alternative bases should be one (not indel)
        if len(ref_base) == 1 and len(alt_base) == 1:
            
            # Reference base should not be ambigous
            if ref_base != 'N':
                
                # Enforce minimum isolate filter call ratio
                if (refRatio >= refRequiredRatio):
                    base_call = 'REF'
                elif refRatio <= altRequiredRatio:
                    base_call = 'ALT'
                else:
                    base_call = 'AMB'

                if base_call == 'ALT':
                    string = isolate_name + '\t'+ isolate_name + '\t'+ str(position) + '\t'+ alt_base + '\t'+ ref_base
                    string_list = [isolate_name,isolate_name,position,alt_base,ref_base]
                    reviewed_SNP_strings.append(string)
                    reviewed_SNP_list.append(string_list)
                    positions_only.append(position)
                    
    print("Number of SNPs after filtering:  ", len(set(positions_only)) )

    reviewed_SNP_df        = pd.DataFrame(reviewed_SNP_list, columns=['isolate_name1','isolate_name2','position','alt_base','ref_base'])
    unique_reviewed_SNP_df = reviewed_SNP_df.drop_duplicates(subset=['position'])

    with open(snp_filter_output, 'w') as outfile:
        for each_string in reviewed_SNP_strings:
            outfile.write(each_string + '\n')


    # Write final SNP list to file
    final_snp_call_list = unique_reviewed_SNP_df.values.tolist()

    with open(final_snp_out, 'w') as final_out:
        for final_snp in final_snp_call_list:
            final_string = str(final_snp[2])
            final_out.write(final_string + '\n')

    return reviewed_SNP_list

#rev_snps = filterSNPs('conf.txt')
#print(len(rev_snps))

