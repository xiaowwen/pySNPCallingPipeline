import numpy as np
import pandas as pd

class parse_vcf(object):

    def __init__(self, vcf_file):

        vcf_file_handle = open(vcf_file)
        header_lines    =  []
        vcf_data_lines  =  []

        for snp_line in vcf_file_handle:
            if not snp_line.startswith('#'):
                data       = snp_line.strip().split()
                chromosome = data[0]
                position   = data[1]
                snp_id     = data[2]
                ref_base   = data[3]
                alt_base   = data[4]
                quality    = data[5]
                vfilter    = data[6]
                info       = data[7]
                vformat    = data[8]
                na00001    = data[9]

                info_split = info.split(";")

                for item in info_split:
                    if item.startswith('DP4'):
                        ref_for, ref_rev, alt_for, alt_rev = item[4:].strip().split(',')

                if info_split[0].startswith('DP'):
                    depth = info_split[0][3:]
                else:
                    depth = 'INDEL'
                vcf_data_lines.append([chromosome,position,ref_base,alt_base,quality,depth,ref_for,ref_rev,alt_for,alt_rev,info,vformat,na00001])

            elif snp_line.startswith('#'):
                header_lines.append(snp_line.strip())

        self.VCF_DataFrame = pd.DataFrame(vcf_data_lines, columns = ['chromosome','position','ref_base','alt_base','quality','depth','ref_for','ref_rev','alt_for','alt_rev','info','vformat','na00001'])
        self.header        = header_lines
