#!/usr/bin/python
import sys
import pandas as pd
import numpy as np
import scipy.stats as st
import math

import seaborn as sns
import matplotlib.pyplot as plt


### Then column standard deviation is not needed?

with open (sys.argv[1], 'r' ) as euk_diso_profile, open(sys.argv[2], 'r') as bacteria_diso_profile, open(sys.argv[3], 'r') as virus_diso_profile:
    euk_diso_profile=pd.read_csv(euk_diso_profile)
    euk_col_means=[]
    euk_col_sd=[]
    euk_col_name=[]
    euk_LI_UI=[]
    for (col_name, col_data) in euk_diso_profile.iteritems():
        df_len=len(col_data)
        
        col_data=[float(x) for x in col_data]
        col_mean=np.nanmean(col_data)
        col_sd=np.nanstd(col_data, ddof=1)
        t_critical=st.t.ppf(q=0.95, df=df_len-1)
        sigma= col_sd/math.sqrt(df_len)
        margin_of_error=t_critical * sigma
        euk_col_means.append(col_mean)
        euk_col_sd.append(col_sd)
        euk_LI_UI.append([col_mean - margin_of_error, col_mean + margin_of_error])
        euk_col_name.append(col_name)

   
    euk_LI=[x[0] for x in euk_LI_UI]
    euk_UI=[x[1] for x in euk_LI_UI]
                            ##############################################################
    bacteria_diso_profile=pd.read_csv(bacteria_diso_profile)
    bac_col_means=[]
    bac_col_sd=[]
    bac_col_name=[]
    bac_LI_UI=[]
    for (col_name, col_data) in bacteria_diso_profile.iteritems():
        df_len=len(col_data)

        col_data=[float(x) for x in col_data]
        col_mean=np.nanmean(col_data)
        col_sd=np.nanstd(col_data, ddof=1)
        t_critical=st.t.ppf(q=0.95, df=df_len-1)
        sigma= col_sd/math.sqrt(df_len)
        margin_of_error=t_critical * sigma
        bac_col_means.append(col_mean)
        bac_col_sd.append(col_sd)
        bac_LI_UI.append([col_mean - margin_of_error, col_mean + margin_of_error])
        bac_col_name.append(col_name)

   
    bac_LI=[x[0] for x in bac_LI_UI]
    bac_UI=[x[1] for x in bac_LI_UI]
                        ###################################################################
    virus_diso_profile=pd.read_csv(virus_diso_profile)
    virus_col_means=[]
    virus_col_sd=[]
    virus_col_name=[]
    virus_LI_UI=[]
    for (col_name, col_data) in virus_diso_profile.iteritems():
        df_len=len(col_data)

        col_data=[float(x) for x in col_data]
        col_mean=np.nanmean(col_data)
        col_sd=np.nanstd(col_data, ddof=1)
        t_critical=st.t.ppf(q=0.95, df=df_len-1)
        sigma= col_sd/math.sqrt(df_len)
        margin_of_error=t_critical * sigma
        virus_col_means.append(col_mean)
        virus_col_sd.append(col_sd)
        virus_LI_UI.append([col_mean - margin_of_error, col_mean + margin_of_error])
        virus_col_name.append(col_name)

    
    virus_LI=[x[0] for x in virus_LI_UI]
    virus_UI=[x[1] for x in virus_LI_UI]           
        
    x=list(range(0,801,1))
    

    fig, ax = plt.subplots(sharey=True,figsize=[20, 10])
    sns.set(font_scale = 3)
    sns.set_style('ticks') 
    
    plt.plot(euk_col_means, color='black')
    plt.plot(euk_LI, color='grey')
    plt.plot(euk_UI, color='grey')
    plt.fill_between(x, y1= euk_UI, y2= euk_LI, color='lightgrey', alpha =0.4)
    
    #x=list(range(200, 401, 1))
    plt.plot(bac_col_means, color='darkgreen')
    plt.plot(bac_LI, color='lightgreen')
    plt.plot(bac_UI, color='lightgreen')
    plt.fill_between(x, y1= bac_UI, y2= bac_LI, color='lightgreen', alpha =0.4)
    
    #x=list(range(300,501,1))
    plt.plot(virus_col_means, color='darkblue')
    plt.plot(x,virus_UI, color='lightblue')
    plt.plot(x,virus_LI, color='lightblue')
    plt.fill_between(x, y1= virus_UI, y2= virus_LI, color='lightblue', alpha =0.4)
    
    
    plt.axvline(400, color='black', linewidth=2, alpha=0.2)
    
    #ax.fill_between(bac_LI, bac_UI, alpha=0.2)
    plt.xticks(list(range(0,801,100)), ['-400','-300', '-200', '-100', '0', '100', '200', '300', '400'], rotation = 90, fontsize=28)
    plt.ylim(0,1)
    ax.set_xlabel("Residues",fontsize=36, fontweight='bold')
    ax.set_ylabel("Short IUPRED2A disorder score",fontsize=36, fontweight='bold')
    ax = plt.gca()
    ax.tick_params(width=3)
    ax.spines['bottom'].set_linewidth(4)
    ax.spines['top'].set_linewidth(4)
    ax.spines['left'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)
    
    #ax.get_figure().savefig('SLiMs_MIDS_profile.png', dpi=600)
