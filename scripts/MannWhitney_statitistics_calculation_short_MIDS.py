from scipy.stats import mannwhitneyu
import sys
import pandas as pd


#calculate Mann-Whitney statistics between different groups
table= sys.argv[1]
with open (table, "r") as analysis_file:
    analysis_table=pd.read_csv(analysis_file)
    
    #create list of the desired groups to test
    
    CLV=analysis_table[(analysis_table.ELMType == 'CLV')].Short_MIDS_per_instance.values.tolist()
    DEG=analysis_table[(analysis_table.ELMType == 'DEG')].Short_MIDS_per_instance.values.tolist()
    DOC=analysis_table[(analysis_table.ELMType == 'DOC')].Short_MIDS_per_instance.values.tolist()
    LIG=analysis_table[(analysis_table.ELMType == 'LIG')].Short_MIDS_per_instance.values.tolist()
    MOD=analysis_table[(analysis_table.ELMType == 'MOD')].Short_MIDS_per_instance.values.tolist()
    TRG=analysis_table[(analysis_table.ELMType == 'TRG')].Short_MIDS_per_instance.values.tolist()
    
    
    
    CLV_euk= analysis_table[(analysis_table.ELMType == 'CLV') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    #print(CLV_euk)

    CLV_bac=analysis_table[(analysis_table.ELMType == 'CLV') & (analysis_table.Taxonomic_group == 'bacteria')].Short_MIDS_per_instance.values.tolist()
    CLV_virus=analysis_table[(analysis_table.ELMType == 'CLV') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    DEG_euk= analysis_table[(analysis_table.ELMType == 'DEG') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    DEG_virus=analysis_table[(analysis_table.ELMType == 'DEG') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    DOC_euk= analysis_table[(analysis_table.ELMType == 'DOC') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    DOC_bac= analysis_table[(analysis_table.ELMType == 'DOC') & (analysis_table.Taxonomic_group == 'bacteria')].Short_MIDS_per_instance.values.tolist()
    DOC_virus=analysis_table[(analysis_table.ELMType == 'DOC') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    LIG_euk= analysis_table[(analysis_table.ELMType == 'LIG') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    LIG_bac= analysis_table[(analysis_table.ELMType == 'LIG') & (analysis_table.Taxonomic_group == 'bacteria')].Short_MIDS_per_instance.values.tolist()
    LIG_virus= analysis_table[(analysis_table.ELMType == 'LIG') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    MOD_euk= analysis_table[(analysis_table.ELMType == 'MOD') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    MOD_bac= analysis_table[(analysis_table.ELMType == 'MOD') & (analysis_table.Taxonomic_group == 'bacteria')].Short_MIDS_per_instance.values.tolist()
    MOD_virus=analysis_table[(analysis_table.ELMType == 'MOD') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    TRG_euk=analysis_table[(analysis_table.ELMType == 'TRG') & (analysis_table.Taxonomic_group == 'eukaryotes')].Short_MIDS_per_instance.values.tolist()
    TRG_bac= analysis_table[(analysis_table.ELMType == 'TRG') & (analysis_table.Taxonomic_group == 'bacteria')].Short_MIDS_per_instance.values.tolist()
    TRG_virus= analysis_table[(analysis_table.ELMType == 'TRG') & (analysis_table.Taxonomic_group == 'viruses')].Short_MIDS_per_instance.values.tolist()
    
    
    #calculate Mann-Whitney statistics
    print('Short disorder IUPRED 2A MIDS comparison by ELM type between taxonomic groups') 
    stat, p_value = mannwhitneyu(CLV_euk, CLV_bac)
    print('CLV_euk vs. CLV_bac' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))

    
    stat, p_value = mannwhitneyu(CLV_euk, CLV_virus)
    print('CLV_euk vs. CLV_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(CLV_bac, CLV_virus)
    print('CLV_bac vs. CLV_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    #####################################NO DEG IN BACTERIA ##################################### OTHERWISE add more analysis to this group
    stat, p_value = mannwhitneyu(DEG_euk, DEG_virus)
    print('DEG_euk vs. DEG_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    

    
    stat, p_value = mannwhitneyu(DOC_euk, DOC_virus)
    print('DOC_euk vs. DOC_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(DOC_euk, DOC_bac)
    print('DOC_euk vs. DOC_bac' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
  
    stat, p_value = mannwhitneyu(DOC_bac, DOC_virus)
    print('DOC_bac vs. DOC_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))



    stat, p_value = mannwhitneyu(LIG_euk, LIG_virus)
    print('LIG_euk vs. LIG_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(LIG_euk, LIG_bac)
    print('LIG_euk vs. LIG_bac' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
  
    stat, p_value = mannwhitneyu(LIG_bac, LIG_virus)
    print('LIG_bac vs. LIG_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    
    stat, p_value = mannwhitneyu(MOD_euk, MOD_virus)
    print('MOD_euk vs. MOD_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(MOD_euk, MOD_bac)
    print('MOD_euk vs. MOD_bac' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
  
    stat, p_value = mannwhitneyu(MOD_bac, MOD_virus)
    print('MOD_bac vs. MOD_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))   
    
    
    stat, p_value = mannwhitneyu(TRG_euk, TRG_virus)
    print('TRG_euk vs. TRG_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(TRG_euk, TRG_bac)
    print('TRG_euk vs. TRG_bac' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))
  
    stat, p_value = mannwhitneyu(TRG_bac, TRG_virus)
    print('TRG_bac vs. TRG_virus' + ' Statistics=%.2f, p=%.6f' % (stat, p_value))   
    
    print('\n Short disorder IUPRED 2A MIDS Virus analysis') 
    ##Virus analysis
    stat, p_value = mannwhitneyu(CLV_virus, DEG_virus)
    print('CLV_virus vs. DEG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(CLV_virus, DOC_virus)
    print('CLV_virus vs. DOC_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_virus, LIG_virus)
    print('CLV_virus vs. LIG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))      
    
    
    stat, p_value = mannwhitneyu(CLV_virus, MOD_virus)
    print('CLV_virus vs. MOD_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_virus, TRG_virus)
    print('CLV_virus vs. TRG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_virus, DOC_virus)
    print('DEG_virus vs. DOC_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(DEG_virus, LIG_virus)
    print('DEG_virus vs. LIG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_virus, MOD_virus)
    print('DEG_virus vs. MOD_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_virus, TRG_virus)
    print('DEG_virus vs. TRG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_virus, LIG_virus)
    print('LIG_virus vs. DOC_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_virus, MOD_virus)
    print('DOC_virus vs. MOD_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_virus, TRG_virus)
    print('DOC_virus vs. TRG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(LIG_virus, MOD_virus)
    print('LIG_virus vs. MOD_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(LIG_virus, TRG_virus)
    print('LIG_virus vs. TRG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(MOD_virus, TRG_virus)
    print('MOD_virus vs. TRG_virus' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    print('\n Short disorder IUPRED 2A MIDS Eukaryotes analysis')    
    ##euk analysis
    stat, p_value = mannwhitneyu(CLV_euk, DEG_euk)
    print('CLV_euk vs. DEG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(CLV_euk, DOC_euk)
    print('CLV_euk vs. DOC_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_euk, LIG_euk)
    print('CLV_euk vs. LIG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))      
    
    
    stat, p_value = mannwhitneyu(CLV_euk, MOD_euk)
    print('CLV_euk vs. MOD_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_euk, TRG_euk)
    print('CLV_euk vs. TRG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_euk, DOC_euk)
    print('DEG_euk vs. DOC_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))
    
    stat, p_value = mannwhitneyu(DEG_euk, LIG_euk)
    print('DEG_euk vs. LIG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_euk, MOD_euk)
    print('DEG_euk vs. MOD_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG_euk, TRG_euk)
    print('DEG_euk vs. TRG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_euk, LIG_euk)
    print('LIG_euk vs. DOC_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_euk, MOD_euk)
    print('DOC_euk vs. MOD_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_euk, TRG_euk)
    print('DOC_euk vs. TRG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(LIG_euk, MOD_euk)
    print('LIG_euk vs. MOD_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(LIG_euk, TRG_euk)
    print('LIG_euk vs. TRG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(MOD_euk, TRG_euk)
    print('MOD_euk vs. TRG_euk' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    print('\n Short disorder IUPRED 2A MIDS Bacteria analysis')        
    ##bac analysis
    
    stat, p_value = mannwhitneyu(CLV_bac, DOC_bac)
    print('CLV_bac vs. DOC_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_bac, LIG_bac)
    print('CLV_bac vs. LIG_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))      
    
    
    stat, p_value = mannwhitneyu(CLV_bac, MOD_bac)
    print('CLV_bac vs. MOD_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV_bac, TRG_bac)
    print('CLV_bac vs. TRG_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    
    stat, p_value = mannwhitneyu(DOC_bac, LIG_bac)
    print('LIG_bac vs. DOC_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_bac, MOD_bac)
    print('DOC_bac vs. MOD_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC_bac, TRG_bac)
    print('DOC_bac vs. TRG_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(LIG_bac, MOD_bac)
    print('LIG_bac vs. MOD_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(LIG_bac, TRG_bac)
    print('LIG_bac vs. TRG_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(MOD_bac, TRG_bac)
    print('MOD_bac vs. TRG_bac' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    


################ All groups by ELM type ###########################
    print('\n Short disorder IUPRED 2A MIDS all analysis')  

    stat, p_value = mannwhitneyu(CLV, DEG)
    print('CLV vs. DEG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV, DOC)
    print('CLV vs. DOC' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV, MOD)
    print('CLV vs. MOD' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(CLV, LIG)
    print('CLV vs. LIG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(CLV, TRG)
    print('CLV vs. TRG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG, DOC)
    print('DEG vs. DOC' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG, LIG)
    print('DEG vs. LIG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG, MOD)
    print('DEG vs. MOD' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DEG, TRG)
    print('DEG vs. TRG' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    
    stat, p_value = mannwhitneyu(DOC, LIG)
    print('DOC vs. LIG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC, MOD)
    print('DOC vs. MOD' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(DOC, TRG)
    print('DOC vs. TRG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(LIG, MOD)
    print('LIG vs. MOD' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(LIG, TRG)
    print('LIG vs. TRG' + ' Statistics=%.2f, p=%f' % (stat, p_value))  
    
    stat, p_value = mannwhitneyu(MOD, TRG)
    print('MOD vs. TRG' + ' Statistics=%.2f, p=%f' % (stat, p_value)) 
    