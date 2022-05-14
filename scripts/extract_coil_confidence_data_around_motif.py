#!/usr/bin/python
import sys
import pandas as pd

RANGE = 100
#This script takes two arguments ex. virus_TP_full_data.csv that has the full data and the column CSV file generated from the extract_coil_confidence_netsurfp ex. virus_netsurfp_coil_conf_column.csv
#The output is 2 CSV files one for the coil confidence of residues before and another for the residues after the motif instance based on the range specified at the begining of the file


ELM_motifs= sys.argv[1]
coil_conf_data=sys.argv[2]

with open (ELM_motifs, 'r') as ELM_motifs, open(coil_conf_data, 'r') as coil_conf_data:
    
    ELM_motifs=pd.read_csv(ELM_motifs).values.tolist()
    coil_conf_data=pd.read_csv(coil_conf_data)
    coil_conf_data=coil_conf_data.groupby(by= 'id', axis=0)

    #print(coil_conf_data)
    
    before_motif=[]
    after_motif=[]
    for key, value in coil_conf_data:
        protein_id= key
        coil_conf_list= value['p[q3_C]'].tolist()
        #print(key)
        #print(coil_conf_list)
      
        for elm in ELM_motifs:
            elm_accession=elm[0]
            elm_type=elm[1]
            elm_identifier=elm[2]
            elm_protein_name=elm[3]
            elm_protein_id= elm[4]
            elm_mean=elm[20]
            #print(elm_protein_id)
            orig_start= int(elm[6])
            orig_end=int(elm[7])
            if protein_id.split('|')[1] == elm_protein_id:
                elm_len = orig_end - orig_start   #Iguess we don't need that line
                
                
                if orig_start - RANGE <= 0:
                    coil_conf = coil_conf_list[0:orig_start-1]
                    #print(orig_start)
                    #print(coil_conf)
                    #print(len(coil_conf))
                    remaining_data_num=RANGE-len(coil_conf)
                    remaining_data=['nan'] * remaining_data_num
                    #print(remaining_data_num)
                    #print(remaining_data)
                    coil_conf_data_upto_RANGE= remaining_data+coil_conf
                    #print(coil_conf_data_upto_RANGE)
                    #print(len(coil_conf_data_upto_RANGE))
                    before_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, coil_conf_data_upto_RANGE])
                   
                else:
                    coil_conf=coil_conf_list[orig_start-(RANGE+1):orig_start-1]
                    #print(orig_start)
                    #print(len(coil_conf))
                    #print(coil_conf)
                    before_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, coil_conf])
                    
                if len(coil_conf_list) - orig_end - RANGE <= 0:
                    
                    coil_conf= coil_conf_list[orig_end:len(coil_conf_list)+1]
                    #print(coil_conf)
                    #print(orig_end)
                    #print(len(coil_conf))
                    remaining_data_num=RANGE-len(coil_conf)
                    remaining_data= ['nan'] * remaining_data_num
                    #print(remaining_data_num)
                    #print(remaining_data)
                    coil_conf_data_upto_RANGE= coil_conf +remaining_data
                    #print(coil_conf_data_upto_RANGE)
                    #print(len(coil_conf_data_upto_RANGE))
                    after_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, coil_conf_data_upto_RANGE])
                else:
                    
                    coil_conf_data_after=coil_conf_list[orig_end:orig_end+RANGE]
                    after_motif.append([elm_accession, elm_type, elm_identifier, elm_protein_name, elm_protein_id, orig_start, orig_end, elm_mean, coil_conf_data_after])
                
                
                
                
                
                
            
    #print(before_motif)
    #print(len(before_motif))  
    before_motif=pd.DataFrame(before_motif)
    before_motif.columns= ['elm_accession', 'elm_type', 'elm_identifier', 'elm_protein_name', 'elm_protein_id', 'orig_start', 'orig_end', 'elm_mean', 'Coil_conf_data_before']
    before_motif.to_csv(str(sys.argv[1]).split('_')[0]+'_coil_conf_before_motif_data_'+str(RANGE)+'.csv', encoding='utf-8', index=False, header= True)
    #print(after_motif)
    #print(len(after_motif))  
    after_motif=pd.DataFrame(after_motif)
    after_motif.columns= ['elm_accession', 'elm_type', 'elm_identifier', 'elm_protein_name', 'elm_protein_id', 'orig_start', 'orig_end', 'elm_mean', 'Coil_conf_data_after']
    after_motif.to_csv(str(sys.argv[1]).split('_')[0]+'_coil_conf_after_motif_data_'+str(RANGE)+'.csv', encoding='utf-8', index=False, header= True)
         
                  
          
                
                
        
                

    