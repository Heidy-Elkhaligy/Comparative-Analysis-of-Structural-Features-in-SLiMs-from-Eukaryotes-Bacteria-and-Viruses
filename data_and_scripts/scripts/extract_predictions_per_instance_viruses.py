#!/usr/bin/python
import sys
import pandas as pd
from statistics import mean

#This script takes the taxonomic group TP instances csv file, netsurfp csv file output, and iupred disorder output
#Output: 1 file merged predictions data / instances in each taxonomic group


with open (sys.argv[1], 'r') as elmfile,open (sys.argv[2],'r') as netsurfp_csv, open (sys.argv[3], 'r') as long_iupred_diso, open(sys.argv[4], 'r') as short_iupred_diso:
    
    #elms=csv.DictReader(elmfile, delimiter=",", quotechar='"')
    elms=pd.read_csv(elmfile)
    column_names=elms.columns.values.tolist()
    elms=elms.values.tolist()
    
    netsurfp_csv=pd.read_csv(netsurfp_csv)
    #print(netsurfp_csv)
    netsurfp_id_group=netsurfp_csv.groupby(by='id', axis=0)
    netsurfp_data=[]

    for key , value in netsurfp_id_group:

        protein_id= key.split('|')[1] # ('|')
        protein_name= key.split('|')[2]
        protein_start=key.split('|')[-1].split('_')[0]
        protein_end=key.split('|')[-1].split('_')[1]
        coil_conf= list(value['p[q3_C]'])
        accessibility=list(value['rsa'])
        ss_pred=list(value['q3'])
        aa=list(value['seq'])
        pos=list(value['n'])
        netsurfp_data.append([protein_id,protein_name,protein_start, protein_end, aa,pos,coil_conf, accessibility, ss_pred])
    #print(netsurfp_data)  
    

    
    
    for elm in elms:
        for protein in netsurfp_data:
            if protein[0]==elm[4] and elm[6] in list(range(int(protein[2]),int(protein[3]))):
                
                if int(protein[2]) == 1:      
                    #print(protein[0])
                    #print(elm[4])
                   # print(elm[0])
                    #print(protein[2])
                    elm_seq_start= elm[6]-1
                    elm_seq_end= elm[7]

                    protein_seq=protein[4][elm_seq_start: elm_seq_end]

                    protein_seq=''.join(protein_seq)
                    elm.append(protein_seq) # append the sequence
                    
                    #accessibility
                    acc=protein[7][elm_seq_start: elm_seq_end] # extract accessibility score
                    acc_bin=['E' if float(x) > 0.25 else 'B' for x in acc] # convert accessibility score to E exposed or B buried
                    acc_bin=''.join(acc_bin)
                    elm.append(acc_bin) # append the accesibility binary values (E/B)
    
                    elm.append(acc_bin.count("E") / len(acc_bin)) # append the accessibility fraction
                    
                    #secondary structure
                    q3=protein[8][elm_seq_start: elm_seq_end]
                    q3_state=''.join(q3)
                    elm.append(q3_state) #append the secondary structure prediction
    
                    elm.append(q3.count("C")/len(q3)) # append the secondary structure fraction
    
                    coil_conf=protein[6][elm_seq_start: elm_seq_end]
                    elm.append(coil_conf)
                    MCCS=mean(coil_conf)
                    elm.append(MCCS) #append MCCS/ instance
                else:
                    #print(protein[0])
                    #print(elm[4])
                    #print(elm[0])
                    #print(protein[2])
                    elm_seq_start= elm[6]-int(protein[2])
                    elm_seq_end= elm[7]-int(protein[2])+1
                    protein_seq=protein[4][elm_seq_start: elm_seq_end]
                    protein_seq=''.join(protein_seq)
                    elm.append(protein_seq) # append the sequence

                    #accessibility
                    acc=protein[7][elm_seq_start: elm_seq_end] # extract accessibility score
                    acc_bin=['E' if float(x) > 0.25 else 'B' for x in acc] # convert accessibility score to E exposed or B buried
                    acc_bin=''.join(acc_bin)

                    elm.append(acc_bin) # append the accesibility binary values (E/B)
    
                    elm.append(acc_bin.count("E") / len(acc_bin)) # append the accessibility fraction
    
                    
                    #secondary structure
                    q3=protein[8][elm_seq_start: elm_seq_end]
                    q3_state=''.join(q3)
                    elm.append(q3_state) #append the secondary structure prediction
    
                    elm.append(q3.count("C")/len(q3)) # append the secondary structure fraction
    
                    
                    coil_conf=protein[6][elm_seq_start: elm_seq_end]
                    elm.append(coil_conf)
                    MCCS=mean(coil_conf)
                    elm.append(MCCS) #append MCCS/ instance
                    
    #Adding disorder analysis
    long_iupred_diso=[line.split() for line in long_iupred_diso]
    long_iupred_diso=[list(item[0].split(',')) for item in long_iupred_diso]    
    
    
    for elm in elms:        
        for value in long_iupred_diso:
            adjusted_start=value[0].split('|')[-1].split('_')[0]

            adjusted_end=value[0].split('|')[-1].split('_')[1]
            if value[0].split('|')[1] == elm[4] and elm[6] in list(range(int(adjusted_start),int(adjusted_end))): 
                if int(adjusted_start) == 1:

                    elm_seq_start= elm[6]
                    elm_seq_end= elm[7]+1
                    diso_values=value[elm_seq_start:elm_seq_end]
                    diso_values=[float(x) for x in diso_values]
    
                    # using 0.4 cutoff
                    diso_04_bin=['O' if x < 0.4 else 'D' for x in diso_values]
                    diso_04_bin=''.join(diso_04_bin)
                    elm.append(diso_04_bin) # append the disorder binary (D/O) based on 0.4 cutoff
    
                    elm.append(diso_04_bin.count('D')/len(diso_04_bin)) # append 0.4 disorder fraction
    
                    
                    #using 0.5 cutoff
                    diso_05_bin=['O' if x < 0.5 else 'D' for x in diso_values]
                    diso_05_bin=''.join(diso_05_bin)
                    elm.append(diso_05_bin) # append the disorder binary (D/O) based on 0.5 cutoff
    
                    elm.append(diso_05_bin.count('D')/len(diso_05_bin)) # append 0.5 disorder fraction
    
                    #MIDS calcuations
                    MIDS=mean(diso_values)
                    elm.append(diso_values)
                    elm.append(MIDS) #append MIDS/instance
                else:
                    #print(value[0].split('|')[1])
                    #print(elm[4])
                    #print(elm[0])
                    #print(adjusted_start)
                    elm_seq_start= elm[6]-int(adjusted_start)+1

                    elm_seq_end= elm[7]-int(adjusted_start)+2

                    diso_values=value[elm_seq_start:elm_seq_end]

                    diso_values=[float(x) for x in diso_values]
    
                    # using 0.4 cutoff
                    diso_04_bin=['O' if x < 0.4 else 'D' for x in diso_values]
                    diso_04_bin=''.join(diso_04_bin)
                    elm.append(diso_04_bin) # append the disorder binary (D/O) based on 0.4 cutoff
    
                    elm.append(diso_04_bin.count('D')/len(diso_04_bin)) # append 0.4 disorder fraction
    
                    
                    #using 0.5 cutoff
                    diso_05_bin=['O' if x < 0.5 else 'D' for x in diso_values]
                    diso_05_bin=''.join(diso_05_bin)
                    elm.append(diso_05_bin) # append the disorder binary (D/O) based on 0.5 cutoff
    
                    elm.append(diso_05_bin.count('D')/len(diso_05_bin)) # append 0.5 disorder fraction
    
                    #MIDS calcuations
                    MIDS=mean(diso_values)
                    elm.append(diso_values)
                    elm.append(MIDS) #append MIDS/instance
                    
    short_iupred_diso=[line.split() for line in short_iupred_diso]
    short_iupred_diso=[list(item[0].split(',')) for item in short_iupred_diso]                    
                    
    for elm in elms:        
        for value in short_iupred_diso:
            adjusted_start=value[0].split('|')[-1].split('_')[0]

            adjusted_end=value[0].split('|')[-1].split('_')[1]
            if value[0].split('|')[1] == elm[4] and elm[6] in list(range(int(adjusted_start),int(adjusted_end))): 
                if int(adjusted_start) == 1:

                    elm_seq_start= elm[6]
                    elm_seq_end= elm[7]+1
                    diso_values=value[elm_seq_start:elm_seq_end]
                    diso_values=[float(x) for x in diso_values]
    
                    # using 0.4 cutoff
                    diso_04_bin=['D' if x > 0.4 else 'O' for x in diso_values]
                    diso_04_bin=''.join(diso_04_bin)
                    elm.append(diso_04_bin) # append the disorder binary (D/O) based on 0.4 cutoff
    
                    elm.append(diso_04_bin.count('D')/len(diso_04_bin)) # append 0.4 disorder fraction
    
                    
                    #using 0.5 cutoff
                    diso_05_bin=['D' if x > 0.5 else 'O' for x in diso_values]
                    diso_05_bin=''.join(diso_05_bin)
                    elm.append(diso_05_bin) # append the disorder binary (D/O) based on 0.5 cutoff
    
                    elm.append(diso_05_bin.count('D')/len(diso_05_bin)) # append 0.5 disorder fraction
    
                    #MIDS calcuations
                    MIDS=mean(diso_values)
                    elm.append(diso_values)
                    elm.append(MIDS) #append MIDS/instance
                else:
                    #print(value[0].split('|')[1])
                    #print(elm[4])
                    #print(elm[0])
                    #print(adjusted_start)
                    elm_seq_start= elm[6]-int(adjusted_start)+1

                    elm_seq_end= elm[7]-int(adjusted_start)+2

                    diso_values=value[elm_seq_start:elm_seq_end]

                    diso_values=[float(x) for x in diso_values]
    
                    # using 0.4 cutoff
                    diso_04_bin=['D' if x > 0.4 else 'O' for x in diso_values]
                    diso_04_bin=''.join(diso_04_bin)
                    elm.append(diso_04_bin) # append the disorder binary (D/O) based on 0.4 cutoff
    
                    elm.append(diso_04_bin.count('D')/len(diso_04_bin)) # append 0.4 disorder fraction
    
                    
                    #using 0.5 cutoff
                    diso_05_bin=['D' if x > 0.5 else 'O' for x in diso_values]
                    diso_05_bin=''.join(diso_05_bin)
                    elm.append(diso_05_bin) # append the disorder binary (D/O) based on 0.5 cutoff
    
                    elm.append(diso_05_bin.count('D')/len(diso_05_bin)) # append 0.5 disorder fraction
    
                    #MIDS calcuations
                    MIDS=mean(diso_values)
                    elm.append(diso_values)
                    elm.append(MIDS) #append MIDS/instance

    

    #print(elms)
    #print(column_names)
    new_names=['SLiM_seq', 'Accessibility_bin', 'Accessibility_fraction', 'SS_prediction', 'Coil_fraction', 'Coil_confidence',\
               'MCCS_per_instance' , 'Long_IUPRED_0.4_bin', 'Long_IUPRED_0.4_fraction', 'Long_IUPRED_0.5_bin', 'Long_IUPRED_0.5_fraction', 'Long_IUPRED_scores','Long_MIDS_per_instance',\
                   'Short_IUPRED_0.4_bin', 'Short_IUPRED_0.4_fraction', 'Short_IUPRED_0.5_bin', 'Short_IUPRED_0.5_fraction', 'Short_IUPRED_scores','Short_MIDS_per_instance' ] 
    merged_data=pd.DataFrame(elms)
   
    merged_data.columns=column_names+new_names
    #drop the rows that contains the SLiMs that doesn't have annotation for their polyproteins
    merged_data=merged_data[merged_data['Primary_Acc'].str.contains('Q8Q7J5-1|Q6PMW3|P03300')==False]

    #print(merged_data)
    
    merged_data.to_csv(str(sys.argv[1]).split('.')[0]+'_full_predictions.csv', encoding='utf-8', index=False, header= True)
    #merged_data= netsurfp_data.apply(lambda x: x.explode() if x.name in ['Coil_conf', 'aa_seq', 'pos'] else x)
    #merged_data.to_csv(str(sys.argv[1]).split('_')[0]+'_netsurfp_coil_conf_column.csv', encoding='utf-8', index=False, header= True)
    
    
   
