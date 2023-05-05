#!/usr/bin/env python
# coding: utf-8
#function
import gzip
import regex
import argparse


def main():
    parser = argparse.ArgumentParser(description='A function to extract gRNA and clonal barcode from merged fastq gz file')
    parser.add_argument("--a", required=True, help="This is the input fastq gz file")
    parser.add_argument("--o", required=True, help="This is the dir of output file")
    args = parser.parse_args()
    fastqgz_input_address = args.a
    output_dir = args.o
    gRNA_output_address = output_dir + '/' + 'gRNA.bartender'
    clonal_barcode_output_address = output_dir + '/' + 'clonalbarcode.bartender'
    
    file_a = open(gRNA_output_address,'wt')
    file_b = open(clonal_barcode_output_address,'wt')
    with gzip.open(fastqgz_input_address,'rt') as handler:
        temp_readID = handler.readline().rstrip() # read ID
        temp_sequence = handler.readline().rstrip()
        handler.readline() # skip two lines
        handler.readline()
        while temp_readID:
            # This is the regular expression pattern
            # 16 bp for barcode
            # 16-20 for sgRNA, some of the control sgRNA are shorter
            temp_pattern = regex.compile('(TAGTT){e<2}' + '(.{16})' + 'TATGG'+'(.{16,21})' + 'GTT(TAAGA){e<2}')

            temp_search_result = temp_pattern.search(temp_sequence)
            if temp_search_result:
                temp_gRNA = temp_search_result.group(3)
                temp_clonal_barcode = temp_search_result.group(2)
                temp_string = '{},{}\n'.format(temp_gRNA,temp_readID) # output to bartender format
                file_a.write(temp_string)
                temp_string = '{},{}\n'.format(temp_clonal_barcode,temp_readID)
                file_b.write(temp_string)
            temp_readID = handler.readline().rstrip() # read ID
            temp_sequence = handler.readline().rstrip()
            handler.readline() # skip two lines
            handler.readline()
    file_a.close()
    file_b.close()
    
if __name__ == "__main__":
    main()  


