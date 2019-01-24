import pandas as pd
from   read_config import read_config
from   parse_vcf import parse_vcf
import glob
import utilities

def write_vcf(config_file):

    # Read configuration file and obtain neccesary file paths, handles and parameters
    config_reader     = read_config(config_file)
    config_dict       = config_reader.get_config_contents()

    reference          = config_dict['REF_FILE']
    sequences          = config_dict['INPUT_DIR']
    alignment_dir      = config_dict['ALIGNMENT_DIR']
    snp_checker_out    = config_dict['SNP_CHECKER_OUTPUT']
    snp_align_fasta    = config_dict['CREATE_ALIGN_OUTPUT']
    final_snp_out      = config_dict['CREATE_ALIGN_LIST_OUTPUT']

    final_SNPs = open(final_snp_out)
    snpchecker = open(snp_checker_out)

    output_vcf = alignment_dir + final_snp_out.split('/')[-1].split('.')[0] + '.vcf'

    # Read Information from the output of SNPChecker

    checked_snp_lines  = pd.read_csv(snp_checker_out)
    unique_checked_snp = checked_snp_lines.drop_duplicates(subset=['position']) 


    # Read in each SNP from final output, use SNPChecker info above to write out a VCF file

    with open(output_vcf, 'w') as outfile:
        for snp_line in final_SNPs:

            data        = unique_checked_snp[unique_checked_snp.position == int(snp_line)]
            chromosome  = data.isolate.values.tolist()[0].split('/')[-1]
            SNP_pos     = data.position.values.tolist()[0]
            SNP_info    = data.snp_info.values.tolist()[0] 
            ID          = '.'
            REF         = data.ref_base.values.tolist()[0] 
            ALT         = data.alt_base.values.tolist()[0]
            QUAL        = data.quality.values.tolist()[0]
            FILTER      = '.'
            VFORMAT     = data.vformat.values.tolist()[0]
            NA00001     = data.na00001.values.tolist()[0]

            string = chromosome+ '\t'+ str(SNP_pos)+ '\t'+ ID+ '\t'+ REF+ '\t'+ ALT+ '\t'+ str(QUAL)+ '\t'+ FILTER+ '\t'+  SNP_info + '\t'+ VFORMAT + '\t'+ NA00001
            
            outfile.write(string + '\n')

#write_vcf('conf.intra.txt')

