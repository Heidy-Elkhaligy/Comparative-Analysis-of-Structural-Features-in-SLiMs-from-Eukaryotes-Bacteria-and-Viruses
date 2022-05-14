#!/usr/bin/python
import csv
from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import pandas as pd


# This file takes two arguments the ELM compiled fasta file and the csv file of all ELMs that includes their taxonomic classification
# The output will be the complete fasta for each ELM category based on their taxonomic classification

seqs=[]
with open (sys.argv[1], 'r') as elm_fasta, open(sys.argv[2], 'r') as elm_csv:
     for title, seq in SimpleFastaParser(elm_fasta):
         #print(title)
         #print(seq)
         seqs.append([title,seq.strip()])  
         
    #print(seqs)
     
     inst_df=pd.read_csv(elm_csv)
     #print(inst_df)
     
     bacteria_data=inst_df[(inst_df['Taxonomic_group'] == 'Bacteria') & (inst_df['InstanceLogic'] == 'true positive')]
     bacteria_accessions= list(bacteria_data.Primary_Acc)
     #print(bacteria_data)
     #print(bacteria_accessions)
     
     virus_data=inst_df[((inst_df['Taxonomic_group'] == 'Viruses_DNA') | (inst_df['Taxonomic_group'] == 'Viruses_RNA')) & (inst_df['InstanceLogic'] == 'true positive')]
     virus_accessions= list(virus_data.Primary_Acc)
     
     eukaryote_data=inst_df[((inst_df['Taxonomic_group'] == 'Eukaryote_Vertebrates') | (inst_df['Taxonomic_group'] == 'Eukaryote_Arthropoda') | (inst_df['Taxonomic_group'] == 'Eukaryote_Fungi')| (inst_df['Taxonomic_group'] == 'Eukaryote_mix')| (inst_df['Taxonomic_group'] == 'Eukaryote_Viridiplantae') | (inst_df['Taxonomic_group'] == 'Eukaryote_Nematode') | (inst_df['Taxonomic_group'] == 'Eukaryote_Amoebozoa')) & (inst_df['InstanceLogic'] == 'true positive')]
     eukaryote_accessions= list(eukaryote_data.Primary_Acc)
     
     
     with open ('bacteria_complete.fasta', 'w') as bact_outp:
         for acc in bacteria_accessions:
            for seq in seqs:
                if acc in seq[0]:
                   bact_outp.write('>'+seq[0]+'\n'+seq[1]+'\n')
            
     with open ('viruses_complete.fasta', 'w') as virus_outp:
        for acc in virus_accessions:
           for seq in seqs:
               if acc in seq[0]:
                  virus_outp.write('>'+seq[0]+'\n'+seq[1]+'\n')   
     
     
     with open ('eukaryotes_complete.fasta', 'w') as eukaryote_outp:
         for acc in eukaryote_accessions:
            for seq in seqs:
                if acc in seq[0]:
                   eukaryote_outp.write('>'+seq[0]+'\n'+seq[1]+'\n')   