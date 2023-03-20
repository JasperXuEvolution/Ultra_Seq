#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import argparse

print("This is a scrip to combine bartender output to generate a dataframe")

def merge_bartender_output(input_barcode_address,input_cluster_address):
    temp_df1 = pd.read_csv(input_barcode_address)
    temp_df2 = pd.read_csv(input_cluster_address)
    temp_merge=pd.merge(temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
         left_index=False, right_index=False, sort=True, copy=True, indicator=False,
         validate=None)
    return(temp_merge)

def unique_read_to_cluster_dic(input_df): # the input is a merged df from two barcode.csv and cluster.csv
    temp_UniqueRead_To_ClusterSeq_dic = {}
    for x,y in zip(input_df['Unique.reads'].to_list(),input_df['Center'].to_list()):
        temp_UniqueRead_To_ClusterSeq_dic[x] = y
    return(temp_UniqueRead_To_ClusterSeq_dic)

def merge_two_barcode(input_bartender_address1,input_bartender_address2,
                      input_dic1,input_dic2,barcode1_name,barcode2_name):
    # input_bartender1 and input_bartender2 should have same sequence id order. The read should be matched.
    temp_dic = {} # key is a (bacode1_cluster_seq, barcode2_cluster_seq), value is the count
    with open(input_bartender_address1, 'r') as handler1, open(input_bartender_address2, 'r') as handler2:
        temp1 = handler1.readline().strip()            
        temp2 = handler2.readline().strip()
        while bool(temp1)&bool(temp2):
            temp_sampleID = temp1.split(',')[1].split('|')[1]
            temp_barcode1_cluster = input_dic1.get(temp1.split(',')[0])
            temp_barcode2_cluster = input_dic2.get(temp2.split(',')[0])
            temp_tuple = tuple([temp_barcode1_cluster,temp_barcode2_cluster,temp_sampleID])
            if temp_tuple in temp_dic.keys():
                temp_dic[temp_tuple]+=1
            else:
                temp_dic[temp_tuple]=1
            temp1 = handler1.readline().strip()            
            temp2 = handler2.readline().strip()
    temp_output_df = pd.DataFrame({barcode1_name: [x[0] for x in temp_dic.keys()],
                                   barcode2_name: [x[1] for x in temp_dic.keys()],
                                   'Sample_ID':[x[2] for x in temp_dic.keys()],
                                  'Frequency': list(temp_dic.values())
                                  })
    return(temp_output_df)

def main():
    parser = argparse.ArgumentParser(description='A function to readout bartender output and combine them to df')
    parser.add_argument("--a1", required=True, help="This is the input file of bacode.csv for barcode1")
    parser.add_argument("--a2", required=True, help="This is the input file of cluster.csv for barcode1")
    parser.add_argument("--a3", required=True, help="This is the input file with barcode1 sequence to sequencing id information")
    parser.add_argument("--a4", required=True, help="This is the name of barcode1")
    parser.add_argument("--b1", required=True, help="This is the input file of bacode.csv for barcode2")
    parser.add_argument("--b2", required=True, help="This is the input file of cluster.csv for barcode2")
    parser.add_argument("--b3", required=True, help="This is the input file with barcode2 sequence to sequencing id information")
    parser.add_argument("--b4", required=True, help="This is the name of barcode2")
    parser.add_argument("--o", required=True, help="This is the name of output file")
    args = parser.parse_args()
    
    barcode_dic1 = unique_read_to_cluster_dic(merge_bartender_output(args.a1,args.a2))
    barcode_dic2 = unique_read_to_cluster_dic(merge_bartender_output(args.b1,args.b2))
    
    temp_df = merge_two_barcode(args.a3,args.b3,barcode_dic1,barcode_dic2,args.a4,args.b4)
    temp_df.to_csv(args.o)
    
if __name__ == "__main__":
    main()  


