# This folder contains the necessary files to redo the ELM analysis done in the Comparative analysis of structural features in SLiMs from eukaryotes, bacteria, and viruses

## The folder contains 3 main subfolders:
---------------------------------------
1. scripts folder that contain all the needed scripts to run the analysis.
2. iupred_disorder folder that contains the results of IUPRED 2A runining on the sequence data for each taxonomic group.
3. netsurfp_output folder that contains the compiled csv output file from running NetSurfP locally. 

## It contains 2 main files:
-------------------------
The ELM instances and FASTA files were downloaded on October 25th, 2021

***5 Jupyter notebooks*** are included and are numbered in the order they should be run.
- Please don't run the whole notebook at once as some steps might need you to do something outside the notebook, and carefully read and follow the instructions.

- Please make sure that you are in the main directory that contain the Jupyter notebooks before you start running the cells, and run the cells in the listed order as each step depends on the results from the previous cell. 


***Note*** The eukaryotes disorder (eukaryotes_complete_filtered_long.result and eukaryotes_complete_filtered_short.result) in the iupred_disorder folder, and eukaryotes netsurfp output (eukaryotes_complete_HHsuite_netsurfp.csv) in netsurfp_output folder, as the one used for the analysis in reference, are not added to this repository. If you want to get access to this data, you can run the analysis locally and follow the instructions in the corresponding README.md file in each prediction folder or you can email the authors of the paper to send you the files to have the complete files before running the analysis. 

