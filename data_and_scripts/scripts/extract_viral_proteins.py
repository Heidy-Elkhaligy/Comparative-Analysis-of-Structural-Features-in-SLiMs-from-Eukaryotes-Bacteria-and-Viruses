#!/usr/bin/python
import sys
import requests
from bs4 import BeautifulSoup as soup
import re
import pandas as pd


with open(sys.argv[1], 'r') as virus_elms, open(sys.argv[2], 'r') as virus_complete_fa:
    virus_table=pd.read_csv(virus_elms)
    virus_table=virus_table[['Primary_Acc','Start', 'End']].values.tolist() 
    
    #This sequence search will help in separating each sequence in the file
    seq_search = re.compile("(>.*?)\n([-\w\n]*)", re.DOTALL) 
    
    #  Reading and saving the input file
    virus_complete_fa=virus_complete_fa.read()
    
    seq_headers=[] #List that will contain sequence headers ordered respectively
    all_seq_list=[]  # List that will contain sequences ordered respectively
    
    for i in re.findall(seq_search, virus_complete_fa): # searching for each individual sequence in the file

        seq_headers.append(i[0]) # extracting the header and adding it to header list
        all_seq_list.append(i[1].replace('\n', '')) # extracting the sequence and adding it to sequences lists

    #get how many proteins if it is a virus polyprotein and their corresponding start, end residues and name
    all_chains=[]   
    no_chain_proteins=[]
    for data in virus_table:
        uniprot_id=data[0]
        if uniprot_id == 'P03076':
            uniprot_id =='P0DOJ8'
            start=data[1]
            end=data[2]
            #print(start)
            #print(type(start))
            
            myurl= 'https://www.uniprot.org/uniprot/'+ str(uniprot_id)+'.gff'
                            
            res= requests.get(myurl)
            html_page= res.content
            
            #parse the uniprot page
            pars_html=soup(html_page,'html.parser')
        
            #extractig text from uniprot page
            text = soup.find_all(pars_html, text=True)
            uniprot_page= ''.join(text)
            for line in uniprot_page.split('\n'):
                line_lst=line.split('\t')
                if 'Chain' in line and start in range(int(line_lst[3]),int(line_lst[4])+1) and end in range(int(line_lst[3]),int(line_lst[4])+1):
                    if 'polyprotein' in line_lst[8].split(';')[1] or'Polyprotein' in line_lst[8].split(';')[1] :
                        pass
                    else:
                        all_chains.append([uniprot_id, int(line_lst[3]), int(line_lst[4]), line_lst[8].split(';')[1]])

            
        else:
            start=data[1]
            end=data[2]
            #print(start)
            #print(type(start))
            
            myurl= 'https://www.uniprot.org/uniprot/'+ str(uniprot_id)+'.gff'
                            
            res= requests.get(myurl)
            html_page= res.content
            
            #parse the uniprot page
            pars_html=soup(html_page,'html.parser')
        
            #extractig text from uniprot page
            text = soup.find_all(pars_html, text=True)
            uniprot_page= ''.join(text)
            
            
            if 'Chain' not in uniprot_page:
                #print(uniprot_page)
                for line in uniprot_page.split('\n'):
                    line_lst=line.split('\t')
                    no_chain_proteins.append(line_lst[0])
            else:         
                for line in uniprot_page.split('\n'):
                    line_lst=line.split('\t')
                    if 'Chain' in line and start in range(int(line_lst[3]),int(line_lst[4])+1) and end in range(int(line_lst[3]),int(line_lst[4])+1):
                        if 'polyprotein' in line_lst[8].split(';')[1] or'Polyprotein' in line_lst[8].split(';')[1] :
                            pass
                        else:
                            all_chains.append([line_lst[0], int(line_lst[3]), int(line_lst[4]), line_lst[8].split(';')[1]])
                    
        #print(all_chains)
        #print(len(all_chains))
        #print(no_chain_proteins)
        #print(len(no_chain_proteins))
        
        seq_dict={}
        for protein_data in all_chains:
            cut_start=protein_data[1]
            cut_end=protein_data[2]
            for name, seq in zip(seq_headers,all_seq_list):
                if protein_data[0] == name.split('|')[1]:
                    cut_seq=seq[cut_start-1:cut_end]
                    seq_dict[name+str('|{}_{}').format(cut_start,cut_end)]=cut_seq
        for protein in no_chain_proteins:
            for name,seq in zip(seq_headers,all_seq_list):
                #print(protein)
                #print(name.split('|')[1])
                if protein == name.split('|')[1]:
                    seq_dict[name+str('|{}_{}').format(1,len(seq))]=seq
                    
    
              
    #print(len(seq_dict))
    with open('viruses_cut_polyproteins.fasta', 'w') as output_file:
        for key, value in seq_dict.items():
            output_file.write(f"{key}\n{value}\n")

