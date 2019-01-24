import pandas as pd
import numpy as np
from   read_config import read_config
from   parse_vcf import parse_vcf
import glob
import utilities

def getHQSNPs(config_file):

    # Read configuration file
    config_reader = read_config(config_file)
    config_dict   = config_reader.get_config_contents()
    
    # Obtain neccesary file paths, handles and parameters
    reference              = config_dict['REF_FILE']
    sequences              = config_dict['INPUT_DIR']
    alignment_dir          = config_dict['ALIGNMENT_DIR']
    hq_snp_quality         = float(config_dict['HQ_SNP_QUALITY'])
    hq_snp_depth           = float(config_dict['HQ_SNP_DEPTH'])
    hq_snp_read_balance    = float(config_dict['HQ_SNP_READ_BALANCE'])
    hq_snp_alt2ref_balance = float(config_dict['HQ_SNP_REF_TO_ALT_BALANCE'])
    hq_snp_cluster_size    = float(config_dict['HQ_SNP_CLUSTER_SIZE'])
    hq_snp_dist_to_contig  = float(config_dict['HQ_SNP_DIST_TO_CONTIG_END'])
    hq_snp_list_out        = config_dict['HQ_SNP_LIST_OUTPUT']
    
    isolates = utilities.isolate_folders(sequences)
    #print('Working on isolates: ', isolates)

    combined_HQ_SNPs = []
    for i in range(len(isolates)):

        combined_HQ_SNPs_for_isolate = []

        # Handles for VCF files for the three alignment methods
        bwa_vcf             =  alignment_dir + '/' + isolates[i] + '/'  + isolates[i] + "_sorted_bwa.vcf"
        lastalign_vcf       =  alignment_dir + '/' + isolates[i] + '/'  + isolates[i] +"_sorted_lastalign.vcf"
        novoalign_vcf       =  alignment_dir + '/' + isolates[i] + '/'  + isolates[i] +"_sorted_novoalign.vcf"

        consensus_SNPs_for_isolate = alignment_dir + '/' + isolates[i] + '/'  + isolates[i] +"_consensus.txt"

        vcf_readers   = [parse_vcf(bwa_vcf), parse_vcf(lastalign_vcf), parse_vcf(novoalign_vcf)]
        all_aligners  = [vcf_reader.VCF_DataFrame.values.tolist() for vcf_reader in vcf_readers]

        # Get all indel positions to aid removal of SNPs within #hq_snp_dist_to_contig# base pairs of an indel
        #HQ_SNPs_for_isolate = []
        for aligner_out in all_aligners:
            HQ_SNPs_for_isolate = [] 
            indel_positions = []
            for line in aligner_out:
                if line[5] == 'INDEL':
                    indel_positions.append( int(line[1]) )
            indel_positions.sort()

            # Go through each line in the VCF and perform filtering for HQ SNPs
            for j in range(len(aligner_out)):
                
                line = aligner_out[j]

                quality    = float(line[4])
                refForward = float(line[6])
                refReverse = float(line[7])
                altForward = float(line[8])
                altReverse = float(line[9])

                readDepth  = line[5] # This variable contains the read depth for SNPs and 'INDEL' string for indels
                
                refTotal   = refForward + refReverse
                altTotal   = altForward + altReverse
                #refRatio   = (refTotal)/(altTotal+refTotal)
                refRatio   = refTotal/altTotal

                # Filter variants that are INDELS and below Phred quality
                if (quality >= hq_snp_quality) and (readDepth != 'INDEL'): 
                    
                    # Filter snps with alt forward and reverse supports below required read balance
                    #if altForward >= hq_snp_read_balance and altReverse >= hq_snp_read_balance:
                    if (altForward > hq_snp_read_balance) and (altReverse > hq_snp_read_balance):

                        # Filter SNPs below specified read depth
                        if float(readDepth) >= hq_snp_depth:

                            # Filter reads that are below specified ALT to REF ratio (This is refRatio >= hq_snp_alt2ref_balance, if refRatio = altTotal/altTotal+refTotal ???)
                            if refRatio <= hq_snp_alt2ref_balance:
                                HQ_SNPs_for_isolate.append( int(line[1]) )

            # Remove SNPs within hq_snp_dist_to_contig (150) bases of an INDEL
            # Maybe use set complement to make this faster, circle through all INDELS and add +/- 150 to each, then do complement of that and HQ_SNPs_for_isolate??

            for indel_position in indel_positions:
                for per_isolate_HQ_SNP in HQ_SNPs_for_isolate:
                    if abs(per_isolate_HQ_SNP - indel_position) < hq_snp_dist_to_contig:
                        HQ_SNPs_for_isolate.pop(HQ_SNPs_for_isolate.index(per_isolate_HQ_SNP))
            
            combined_HQ_SNPs_for_isolate.append(HQ_SNPs_for_isolate)

        # Combine HQ SNPs from all three alignment methods and get consensus SNPs

        BWA_aln  = combined_HQ_SNPs_for_isolate[0]
        LAST_aln = combined_HQ_SNPs_for_isolate[1]
        NOVO_aln = combined_HQ_SNPs_for_isolate[2]

        consensus_snps = list( set(BWA_aln) & set(LAST_aln) & set(NOVO_aln) )
        consensus_snps.sort()
        combined_HQ_SNPs.append(consensus_snps)
        np.savetxt(consensus_SNPs_for_isolate, np.array(consensus_snps).astype(int), fmt='%i', delimiter=",")

    SNP_set = set()
    for combined_HQ_SNP in combined_HQ_SNPs:
        SNP_set = SNP_set.union(set(combined_HQ_SNP))
    HQ_SNPs = list(SNP_set)
    HQ_SNPs.sort()

    # Remove SNPs that are within # hq_snp_cluster_size # of each other
    declustered_HQ_SNPs = []
    declustered_HQ_SNPs.append(HQ_SNPs[0])

    for k in range(1,len(HQ_SNPs)):
        check_CLUSTER = ( int(HQ_SNPs[k]) - int(HQ_SNPs[k-1]) ) > hq_snp_cluster_size #>= hq_snp_cluster_size
        if check_CLUSTER:
            declustered_HQ_SNPs.append(HQ_SNPs[k])

    # Remove SNPs that are within hq_snp_dist_to_contig bases of the REFERENCE contigs. This cannot be done with
    # INDELS because this is reference wide

    # First read REFERENCE and get contig ends.
    ref_f = open(reference)
    lines = ref_f.readlines()
    ref_length  = 0
    contig_ends = []
    count_ends  = 0
    for line in lines:
        if line.startswith('>'):
            count_ends += 1
            contig_ends.append(count_ends)
        else:
            ref_length = ref_length + len(line.strip())
    # Add the last reference base 
    contig_ends.append(ref_length)
    ref_f.close()

    # Now remove #the contig_ends# obtained above
    for contig_end in contig_ends:
        for declustered_HQ_SNP in declustered_HQ_SNPs:
            if abs(declustered_HQ_SNP - contig_end) < hq_snp_dist_to_contig:
                declustered_HQ_SNPs.pop(declustered_HQ_SNPs.index(declustered_HQ_SNP))

    SNP_pos_matrix = np.array(declustered_HQ_SNPs)
    np.savetxt(hq_snp_list_out, SNP_pos_matrix.astype(int), fmt='%i', delimiter=",")

    return declustered_HQ_SNPs

#getHQSNPs('conf.intra.txt')
