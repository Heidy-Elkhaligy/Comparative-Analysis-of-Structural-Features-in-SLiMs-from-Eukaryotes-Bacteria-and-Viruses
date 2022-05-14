#!/usr/bin/python
import sys
import pandas as pd
import re

RANGE = sys.argv[3] #chage the range as you want


#This script takes two arguments virus_TP_full_data.csv that has the full data and the matrix file generated from Janelle's disorder prediction script virus_complete_filtered.result_0.5.mtx
#The output is 2 CSV files one for the disorder of residues before and another for the residues after the motif instance based on the range specified at the begining of the file

ELM_motifs= sys.argv[1]
disorder_data=sys.argv[2]

with open (ELM_motifs, 'r') as ELM_motifs, open(disorder_data, 'r') as disorder_data:
    
    ELM_motifs=pd.read_csv(ELM_motifs).values.tolist()
    disorder_data=pd.read_csv(disorder_data,header=None, sep='\r', lineterminator='\n').values.tolist() 
    #disorder_data=disorder_data.values.tolist()
    #print(ELM_motifs)
    #print(disorder_data[0])
    
  
    before_motif=[]
    after_motif=[]

    for index, data in enumerate(disorder_data):
        protein_id= data[0].split('|')[1]
        data_lst=re.split(r'\[|\]|,| ', data[0])[1:]
        #print(data_lst)
        #print(len(data_lst))
        #print(protein_id)
        for elm in ELM_motifs:
            elm_accession=elm[0]
            elm_type=elm[1]
            elm_identifier=elm[2]
            elm_protein_name=elm[3]
            elm_protein_id= elm[4]
            elm_mean=elm[-1]
            #print(elm_protein_id)
            orig_start= int(elm[6])
            orig_end=int(elm[7])
            #print(data)
            if protein_id == elm_protein_id:
                #get the length of the instance
                elm_len = orig_end - orig_start
                #if the start of instance - the range is less than zero then the instance starts before site 401
                if orig_start - RANGE <= 0:

                    diso_data = data_lst[0:orig_start-1]
                    #print(orig_start)
                    #print(diso_data)
                    #print(len(diso_data))

                    remaining_data_num=RANGE-len(diso_data)
                    remaining_data=['nan'] * remaining_data_num
                    #print(remaining_data_num)
                    #print(remaining_data)
                    diso_data_upto_RANGE= remaining_data+diso_data
                    #print(diso_data_upto_RANGE)
                    #print(len(diso_data_upto_RANGE))
                    before_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, diso_data_upto_RANGE])
                   
                else:
                    diso_data=data_lst[orig_start-(RANGE+1):(orig_start-1)]
                    #print(orig_start)
                    #print(len(diso_data))
                    #print(diso_data)
                    before_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, diso_data])
                    
                if len(data_lst) - orig_end - RANGE <= 0:
                    
                    diso_data= data_lst[orig_end:len(data_lst)+1]
                    #print(diso_data)
                    #print(orig_end)
                    #print(len(diso_data))
                   # print(elm_accession)
                    remaining_data_num=RANGE-len(diso_data)
                    remaining_data= ['nan'] * remaining_data_num
                    #print(remaining_data_num)
                    #print(remaining_data)
                    diso_data_upto_RANGE= diso_data + remaining_data
                    #print(diso_data_upto_RANGE)
                    #print(len(diso_data_upto_RANGE))
                    after_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, diso_data_upto_RANGE])
                else:
                    
                    diso_data_after=data_lst[orig_end:orig_end+RANGE]
                    after_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, diso_data_after])
                    #print(diso_data_after)
                    #print(orig_end)
                    #print(len(diso_data_after))
                    #print(elm_accession)
               
                
    #print(before_motif)
    #print(len(before_motif))  
    before_motif=pd.DataFrame(before_motif)
    before_motif.columns= ['elm_accession', 'elm_type', 'elm_identifier', 'elm_protein_name', 'elm_protein_id', 'orig_start', 'orig_end', 'elm_mean', 'diso_data_before']
    before_motif.to_csv(str(sys.argv[1]).split('_')[0]+'_disorder_before_motif_data_'+str(RANGE)+'.csv', encoding='utf-8', index=False, header= True)
    #print(after_motif)
    #print(len(after_motif))  
    after_motif=pd.DataFrame(after_motif)
    after_motif.columns= ['elm_accession', 'elm_type', 'elm_identifier', 'elm_protein_name', 'elm_protein_id', 'orig_start', 'orig_end', 'elm_mean', 'diso_data_after']
    after_motif.to_csv(str(sys.argv[1]).split('_')[0]+'_disorder_after_motif_data_'+str(RANGE)+'.csv', encoding='utf-8', index=False, header= True)
                
                     
                
                
                
                
        
                

    