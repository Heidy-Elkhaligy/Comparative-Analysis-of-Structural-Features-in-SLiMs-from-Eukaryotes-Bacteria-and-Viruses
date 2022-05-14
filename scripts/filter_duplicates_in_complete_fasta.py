import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd

#This script filters the duplicated fasta files and only keeps one
fasta_f= sys.argv[1]

with open( fasta_f, 'r') as complete_fa, open (str(fasta_f).split('.')[0]+'_filtered.fasta', 'w') as outp:
    filtered_fasta=[]
    for title, seq in SimpleFastaParser(complete_fa):
        fasta= [title, seq.strip()]
        #print(fasta)
        if fasta not in filtered_fasta:
            filtered_fasta.append([title,seq.strip()])  
            
    #print(len(filtered_fasta))
    
    for fasta in filtered_fasta:
        outp.write('>'+fasta[0]+'\n'+fasta[1]+'\n')
        
    