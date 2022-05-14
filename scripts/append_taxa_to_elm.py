#!/usr/bin/python

import sys
import csv

with open(sys.argv[1], 'r') as csvfile, open(sys.argv[2], 'r') as tax_data:
    
    elms = csv.DictReader(filter(lambda row: row[0] != "#", csvfile), delimiter="\t", quotechar='"')
    

    Eukaryote_Fungi=[]
    Eukaryote_Nematode=[]
    Eukaryote_Vertebrates=[]
    Eukaryote_Amoebozoa=[]
    Eukaryote_Arthropoda=[]
    Eukaryote_Viridiplantae=[]
    Eukaryote_mix=[]
    
    Bacteria=[]
    
    Viruses_RNA=[]
    Viruses_DNA=[]
    unspecified=[]
    
    #add the taxonomic group to the dataframe
    tax_df = csv.reader(tax_data , quotechar='"')
    for data in tax_df:
        species= data[0].strip().replace('"','')
        
        taxonomy_long=data[1].strip('\t').split(' ')
        #print(species)
        if '2759' and '4751' in taxonomy_long:
            Eukaryote_Fungi.append(species)
            
        elif '2759' and '6231' in taxonomy_long:
            Eukaryote_Nematode.append(species)
            
        elif '2759' and '7742' in taxonomy_long:
            Eukaryote_Vertebrates.append(species)
        
        elif '2759' and '554915' in taxonomy_long:
            Eukaryote_Amoebozoa.append(species)
            
        elif '2759' and '6656' in taxonomy_long:
            Eukaryote_Arthropoda.append(species)
            
        elif '2759' and '33090' in taxonomy_long:
           Eukaryote_Viridiplantae.append(species)
           
        elif '2' in taxonomy_long:
            Bacteria.append(species)
           
        elif '10239' and '2559587' in taxonomy_long:
           Viruses_RNA.append(species)
       
        elif '10239' and '2731341'  in taxonomy_long:
           Viruses_DNA.append(species)
       
        elif '10239' and   '2731342' in taxonomy_long:
           Viruses_DNA.append(species)
           
        elif '10239' and  '2732004'  in taxonomy_long:
            Viruses_DNA.append(species)
           
        elif '10239' and  '2840056' in taxonomy_long:
            Viruses_DNA.append(species)
            
        elif '2759'  in taxonomy_long:
            Eukaryote_mix.append(species)
           
        else:
            unspecified.append([species,'unspecified'])


    #print(unspecified)
    #print(Viruses_RNA)


    ########## Modify the original instances tsv file by adding the taxonomic group to each instance  ################
    with open('elm_instances_with_taxa.csv', 'w', newline='') as csv_file: # change the name of the file
        fieldnames = ['Accession','ELMType','ELMIdentifier','ProteinName','Primary_Acc','Accessions','Start','End','References','Methods','InstanceLogic','PDB','Organism','Taxonomic_group']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        for elm in elms:
           # print(elm['Organism'])
            if elm['Organism'] in Eukaryote_Fungi:
                elm['Taxonomic_group'] = 'Eukaryote_Fungi'
                writer.writerow(elm)
            elif elm['Organism'] in Eukaryote_Nematode:
                elm['Taxonomic_group'] = 'Eukaryote_Nematode'
                writer.writerow(elm)
            elif elm['Organism'] in Eukaryote_Vertebrates:
                elm['Taxonomic_group'] = 'Eukaryote_Vertebrates'
                writer.writerow(elm)
            elif elm['Organism'] in Eukaryote_Amoebozoa:
                elm['Taxonomic_group'] = 'Eukaryote_Amoebozoa'
                writer.writerow(elm)
            elif elm['Organism'] in Eukaryote_Arthropoda:
                elm['Taxonomic_group'] = 'Eukaryote_Arthropoda'
                writer.writerow(elm)
            elif elm['Organism'] in Eukaryote_Viridiplantae:
                elm['Taxonomic_group'] = 'Eukaryote_Viridiplantae'
                writer.writerow(elm)
                
            elif elm['Organism'] in Eukaryote_mix:
                elm['Taxonomic_group'] = 'Eukaryote_mix'
                writer.writerow(elm)
                
            elif elm['Organism'] in Bacteria:
                elm['Taxonomic_group'] = 'Bacteria'
                writer.writerow(elm)
            
            elif elm['Organism'] in Viruses_RNA:
                elm['Taxonomic_group'] = 'Viruses_RNA'
                writer.writerow(elm)
            elif elm['Organism'] in Viruses_DNA:
                elm['Taxonomic_group'] = 'Viruses_DNA'
                writer.writerow(elm)
            else:
                #print(elm['Organism'])
                elm['Taxonomic_group'] = 'unspecified'
                writer.writerow(elm)
            
             