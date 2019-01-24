import pandas as pd
import numpy as np
from   read_config import read_config
from   parse_vcf import parse_vcf
import glob
import utilities

def checkSNPs(config_file):

    # Read configuration file
    config_reader = read_config(config_file)
    config_dict   = config_reader.get_config_contents()
    
    # Obtain neccesary file paths, handles and parameters
    reference              = config_dict['REF_FILE']
    sequences              = config_dict['INPUT_DIR']
    alignment_dir          = config_dict['ALIGNMENT_DIR']
    snp_checker_in         = config_dict['SNP_CHECKER_INPUT']
    snp_filter_quality     = float(config_dict['SNP_FILTER_MIN_GOOD_QUALITY'])
    snp_filter_balance     = float(config_dict['SNP_FILTER_MIN_GOOD_BALANCE'])

    isolates = utilities.isolate_folders(sequences)
    HQ_SNP_list = [int(snp.strip()) for snp in open(snp_checker_in)]
    HQ_SNP_list.sort()

    # Load all isolates, get their VCF objects and put it in a list here so that file IO doesn't happen at every step

    all_isolate_dfs = []
    for i in range(len(isolates)):

        # Read consensus HQ_SNPs from each isolate folder
        bwa_vcf_names  = alignment_dir + '/' + isolates[i] + '/'  + isolates[i] +"_sorted_bwa.vcf"
        vcf_reader     = parse_vcf(bwa_vcf_names)
        isolate_df_    = vcf_reader.VCF_DataFrame
        #all_isolate_dfs.append([isolates[i],isolate_df_])
        all_isolate_dfs.append([bwa_vcf_names,isolate_df_])

    #columns = ['chromosome','position','ref_base','alt_base','quality','depth','ref_for','ref_rev','alt_for','alt_rev','info']
    combined_HQ_SNPs = []
    for HQ_SNP in HQ_SNP_list:

        for isolate_df_info in all_isolate_dfs:

            iso_name   = isolate_df_info[0]
            isolate_df = isolate_df_info[1]
            snp_slice  = isolate_df[isolate_df.position == str(HQ_SNP)]

            if len(snp_slice) != 0:

                line       = snp_slice.values.tolist()[0]
                ref_base   = line[2]
                alt_base   = line[3]
                quality    = float(line[4])
                readDepth  = line[5] # This variable contains the read depth for SNPs and 'INDEL' string for indels
                refForward = float(line[6])
                refReverse = float(line[7])
                altForward = float(line[8])
                altReverse = float(line[9])

                refTotal   = refForward + refReverse
                altTotal   = altForward + altReverse
                refRatio   = refTotal/altTotal
                snp_info   = line[10]
                vformat    = line[11]
                na00001    = line[12]

                # Filter variants, that are INDELS, below Phred quality and required read balance
                if quality >= snp_filter_quality:
                    if readDepth != 'INDEL':
                        if (altForward >= snp_filter_balance) and (altReverse >= snp_filter_balance):
                            combined_HQ_SNPs.append([iso_name, line[1], ref_base, alt_base, quality, readDepth, refForward, refReverse, altForward, altReverse, refRatio, snp_info, vformat, na00001])
                            #combined_HQ_SNPs.append([line[1], ref_base, alt_base, quality, readDepth, refForward, refReverse, altForward, altReverse, refRatio])

    # Write checked SNPs to file

    snp_checker_out    = config_dict['SNP_CHECKER_OUTPUT']
    snp_checker_out_df = pd.DataFrame(combined_HQ_SNPs, columns=['isolate', 'position', 'ref_base', 'alt_base', 'quality', 'readDepth', 'refForward', 'refReverse', 'altForward', 'altReverse', 'refRatio', 'snp_info', 'vformat','na00001'])
    snp_checker_out_df.to_csv(snp_checker_out, index=False)

    #with open(snp_checker_out, 'w') as outfile:
    #    for checked_snp in combined_HQ_SNPs:
    #        string = ''
    #        for item in checked_snp:
    #            string = string + str(item) + ',' + '\t'
    #        outfile.write(string + '\n')

#checkSNPs('conf.txt')

