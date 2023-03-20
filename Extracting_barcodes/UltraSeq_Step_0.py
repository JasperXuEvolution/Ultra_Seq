#!/usr/bin/env python
# coding: utf-8
import gzip
import glob
import pandas as pd
import argparse
def combine_fastq(input_sample_info_address,r1_output_address,r2_output_address):
    # input_sample_info_address contains all the fastq file address
    info_df = pd.read_csv(input_sample_info_address)
     # group based on R1 or R2 read
    info_df_r1 = info_df[info_df['Sample_dir'].str.contains("_1.fq")].sort_values(by='Sample_ID')
    info_df_r2 = info_df[info_df['Sample_dir'].str.contains("_2.fq")].sort_values(by='Sample_ID')
    
    # open output file
    file_r1 = open(r1_output_address,'wt')
    file_r2 = open(r2_output_address,'wt')

    for (index1,x),(index2,y) in zip(info_df_r1.iterrows(),info_df_r2.iterrows()):

        with gzip.open(x['Sample_dir'],'rt') as handler1, gzip.open(y['Sample_dir'],'rt') as handler2:
            #read ID
            temp_readID_r1 = handler1.readline().rstrip() # read ID for R1
            temp_readID_r2 = handler2.readline().rstrip() # read ID for R2
            while temp_readID_r1 and temp_readID_r2:
                #modify id
                temp_id_r1 = temp_readID_r1.split(' ')[0] + '|' + x['Sample_ID'] 
                file_r1.write(temp_id_r1 + '\n')

                temp_id_r2 = temp_readID_r2.split(' ')[0] + '|' + y['Sample_ID']
                file_r2.write(temp_id_r2 + '\n')

                #read and write rest of the file
                file_r1.write(handler1.readline())
                file_r1.write(handler1.readline())
                file_r1.write(handler1.readline())

                file_r2.write(handler2.readline())
                file_r2.write(handler2.readline())
                file_r2.write(handler2.readline())

                temp_readID_r1 = handler1.readline().rstrip() # read ID for R1
                temp_readID_r2 = handler2.readline().rstrip() # read ID for R2         
    file_r1.close()
    file_r2.close()

def main():
    parser = argparse.ArgumentParser(description='A function to combined the fastq files and add Sample_ID information into the description')
    parser.add_argument("--a1", required=True, help="This is the address of sample info")
    parser.add_argument("--o", required=True, help="This the output dir")
    args = parser.parse_args()
    
    sample_info_dir = args.a1
    output_dir = args.o
    # Output address
    combined_fastq_r1_output_address = output_dir + '/Combined_R1.fastq'
    combined_fastq_r2_output_address = output_dir + '/Combined_R2.fastq'
    
    combine_fastq(sample_info_dir,combined_fastq_r1_output_address,combined_fastq_r2_output_address)
    
if __name__ == "__main__":
    main()  
