#!/usr/bin/python

import sys 
import re
import pandas as pd

# This script takes the iupred.result file downloaded from IUPRED 2A 
# and generate a CSV file where each row represent the disorder scores for each protein and the first column contains the protein header

with open(sys.argv[1], 'rt') as long_iupred_result, open(str(sys.argv[1]).split('_')[0]+'_long_diso_scores.csv','w') as long_diso_scores,\
    open(sys.argv[2], 'rt') as short_iupred_result, open(str(sys.argv[2]).split('_')[0]+'_short_diso_scores.csv','w') as short_diso_scores:
      
     long_iupred_result = pd.read_table(long_iupred_result,comment='#', names=['Pos','AA','IUPRED_score', 'ANCHOR_score'], na_values=True, skip_blank_lines=True).fillna(1)
     long_headers_index=long_iupred_result.index[long_iupred_result['Pos'].str.contains('>', na=False)].tolist()
     long_all_seq=[]
     for index in range(len(long_headers_index)):
         try:
             start=long_headers_index[index]
             end=long_headers_index[index+1]
             seq_data=long_iupred_result.loc[range(start,end)]
             #print(seq_data)
             long_all_seq.append(seq_data.values.tolist())
         except IndexError:
             start=long_headers_index[index]
             end=len(long_iupred_result.values.tolist())
             seq_data=long_iupred_result.loc[range(start,end)]
             #print(seq_data.values.tolist())
             long_all_seq.append(seq_data.values.tolist())
     
     #print(all_seq)   
     
     #extract the IUPRED disorder values 
     
     for seq in long_all_seq:
         header=seq[0][0]
         
         seq_diso_values=seq[1:len(seq)]
         diso_04_bin=[]
         diso_05_bin=[]
         diso_value=[]
         for lst in seq_diso_values:
             iupred_diso=lst[2]
             diso_value.append(iupred_diso)
         diso_value=list(map(str,diso_value))
         diso_value=','.join(diso_value)
         long_diso_scores.write(header.split('>')[1]+','+diso_value+'\n')
         
     short_iupred_result = pd.read_table(short_iupred_result,comment='#', names=['Pos','AA','IUPRED_score', 'ANCHOR_score'], na_values=True, skip_blank_lines=True).fillna(1)
     short_headers_index=short_iupred_result.index[short_iupred_result['Pos'].str.contains('>', na=False)].tolist()
     short_all_seq=[]
     for index in range(len(short_headers_index)):
         try:
             start=short_headers_index[index]
             end=short_headers_index[index+1]
             seq_data=short_iupred_result.loc[range(start,end)]
             #print(seq_data)
             short_all_seq.append(seq_data.values.tolist())
         except IndexError:
             start=short_headers_index[index]
             end=len(short_iupred_result.values.tolist())
             seq_data=short_iupred_result.loc[range(start,end)]
             #print(seq_data.values.tolist())
             short_all_seq.append(seq_data.values.tolist())
     
     #print(all_seq)   
     
     #extract the IUPRED disorder values 
     
     for seq in short_all_seq:
         header=seq[0][0]
         
         seq_diso_values=seq[1:len(seq)]
         diso_04_bin=[]
         diso_05_bin=[]
         diso_value=[]
         for lst in seq_diso_values:
             iupred_diso=lst[2]
             diso_value.append(iupred_diso)
         diso_value=list(map(str,diso_value))
         diso_value=','.join(diso_value)
         short_diso_scores.write(header.split('>')[1]+','+diso_value+'\n')