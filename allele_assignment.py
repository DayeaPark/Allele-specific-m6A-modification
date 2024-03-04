
import gzip
import argparse
import re
import pdb

import sys
import os

import numpy as np


script_folder = os.path.dirname(os.path.realpath(__file__))


sys.path.insert(0, os.path.dirname(script_folder) )

import ref_lib
from ref_lib.Vcf import VcfFile

from collections import defaultdict, OrderedDict

import pysam

##################################################################

def get_parameters():
    parser = argparse.ArgumentParser(
                  description = "Reports the count of each nucleotide at each "
                                "nucleotide position. "
                                "Reads with INDELS are discarded.")

    parser.add_argument("--bam",
                        type     = str,
                        required = True,
                        help     = "Bam file coming from transcriptomic alignments")
    
    parser.add_argument("--threshold_ratio",
                        type     = float,
                        required = False,
                        default  = 0.7, 
                        help     = "Minimum ratio threshold for alternative snps.")
    
    parser.add_argument("--threshold_snp",
                        type     = float,
                        required = False,
                        default  = 2, 
                        help     = "Minimum SNP threshold per read.")

    parser.add_argument("--vcf",
                        type     = str,
                        required = True,
                        help     = "vcf file containing the snps")


 
    args = parser.parse_args()
    return args

##############################################################################

def read_vcf(vcf_file):
    vcf_reader = VcfFile(vcf_file)

    contents_dict = OrderedDict()

    for entry in vcf_reader:

        transcript_id = entry.fields["CHROM"]
        position      = int(entry.fields["POS"])
        ref           = entry.fields["REF"]
        alt           = entry.fields["ALT"]

        if not contents_dict.get(transcript_id):
            contents_dict[transcript_id] = OrderedDict()

        contents_dict[transcript_id][position] = { "REF": ref, "ALT": alt,
                                                   "A" : 0, "C" : 0,
                                                   "G" : 0, "T" : 0 , "N" : 0}

    return contents_dict

###########################################################################
### snp position calculation     

def my_mapper(cigar, start_pos, snp_pos, read_length):   
    position = snp_pos - start_pos 
    number_of_matches = 0       
    number_of_insertion = 0
    number_of_deletion = 0
    clip = 0
    
    #print("###start###")
    #print(cigar)
    for index in range(0, len(cigar)):
        if cigar[index][0] in [4,5]:
            if cigar[index][1] >= read_length:
                return -2
            elif cigar[index][1] < read_length:
                clip = clip + cigar[index][1]
                
        elif cigar[index][0] == 0 : # match
            number_of_matches = number_of_matches + cigar[index][1]
            if number_of_matches  + number_of_deletion >= position: 
                if clip + (position - number_of_deletion) + number_of_insertion < read_length:
                    position_in_read  = clip + (position - number_of_deletion) + number_of_insertion
                    return(position_in_read)
                elif number_of_matches  + number_of_deletion >= read_length:
                    return -2
        
        elif cigar[index][0] == 1 : # insertion
            number_of_insertion = number_of_insertion + cigar[index][1]
            
        elif cigar[index][0] == 2 : # deletion
            number_of_deletion = number_of_deletion + cigar[index][1]
            if number_of_matches + number_of_deletion  >= position: 
                #print ("SNP is on a deletion site.")
                return -1
        
        else: 
            #print("Non M,D,I,S,H found!\n")
            return -3
        
    
    return -10
        
          
##############################################################################
        
def count_nucleotides(vcf_contents, bam_file, threshold_ratio, threshold_snp):
    """
    Note that psam coordinates are 0-based and vcf coordinates are 1-based
    So we need to make a conversion
    """
    
    threshold_ref = 1 - threshold_ratio ## default threshold_ratio : 0.7
    
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    reference_file = pysam.AlignmentFile("reference.bam", "wb", template=samfile)
    alternative_file = pysam.AlignmentFile("alternative.bam", "wb", template=samfile)
    undefined_file = pysam.AlignmentFile("undefined.bam", "wb", template=samfile)
    
    for transcript_header, pos_dict in vcf_contents.items():
        
        sorted_positions = sorted(list(pos_dict.keys() ) )
        
        for alignment in samfile.fetch(contig=transcript_header):
            print("###")
            
            read_start_position = alignment.reference_start + 1
            read_end_position = alignment.reference_end + 1

            snp_start_index = np.searchsorted(sorted_positions, read_start_position)
            snp_stop_index  = np.searchsorted(sorted_positions, read_end_position)
                        
            if read_start_position == read_end_position:
                print(alignment)

            if alignment.query_sequence == None:
                print("no read")
                undefined_file.write(alignment)
                continue


            alt_count = 0
            ref_count = 0
            no_match = 0
            unknown = 0
            SNP_no_count = 0
            
            
            for i in range( snp_start_index, snp_stop_index ):
                
                position_in_read = my_mapper(cigar = alignment.cigartuples, 
                                              start_pos = read_start_position, 
                                              snp_pos = sorted_positions[i],
                                              read_length = alignment.query_alignment_length)

                if position_in_read >= 0:
                    nucleotide_in_read = alignment.query_sequence[ position_in_read ]
                
                    if nucleotide_in_read == pos_dict[ sorted_positions[i]]["REF"]:
                        ref_count += 1
                        continue
                    elif nucleotide_in_read == pos_dict[ sorted_positions[i]]["ALT"]:
                        alt_count += 1
                        continue
                    elif nucleotide_in_read != pos_dict[ sorted_positions[i]]["ALT"] and nucleotide_in_read != pos_dict[ sorted_positions[i]]["REF"]:
                        no_match += 1
                        continue
                    else:
                        unknown += 1
                        continue
         
                elif position_in_read < 0:    ## We do not count SNP in deletion site.  
                    SNP_no_count += 1 
            
            print(transcript_header, alt_count, ref_count, unknown)    

            if ref_count + alt_count > threshold_snp: 
                paternal_ratio = alt_count / (ref_count + alt_count) 
                
                if paternal_ratio >= threshold_ratio:
                    alternative_file.write(alignment)
            
                elif paternal_ratio <= threshold_ref:
                    reference_file.write(alignment)
            
                else:
                    undefined_file.write(alignment)
                    
            elif ref_count + alt_count <= threshold_snp:
                undefined_file.write(alignment)      
                
     
    alternative_file.close()    
    reference_file.close()  
    undefined_file.close()  
    
    return vcf_contents

##############################################################################

def main():
    arguments = get_parameters()

    vcf_contents = read_vcf(arguments.vcf)

    nucleotide_counts_dict = count_nucleotides(vcf_contents, arguments.bam, threshold_ratio=arguments.threshold_ratio, threshold_snp=arguments.threshold_snp)
  

if __name__ == "__main__":
    main()
