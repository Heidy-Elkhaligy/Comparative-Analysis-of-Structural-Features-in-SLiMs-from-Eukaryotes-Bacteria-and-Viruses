import pandas as pd


def merge_score_data(data_before, data_after):
    #RANGE=int(RANGE)
    ELM_mean=data_before.elm_mean
    data_before=data_before.iloc[: , -1].str.split(',', expand=True).replace( ["\\[", "\\]", "'"], "", regex=True)
    data_after=data_after.iloc[: , -1].str.split(',', expand=True).replace( ["\\[", "\\]", "'"], "", regex=True)
    merged_data=[data_before,ELM_mean,data_after]
    merged_df=pd.concat(merged_data,axis=1, ignore_index=True)
    merged_df.columns=list(map(str,(list(range(-100,101)))))
    return merged_df

#Eukaryotes Disorder profile for SLiMs and 400 residues before and after
eukaryotes_diso_before=pd.read_csv('viruses_long_disorder_before_motif_data_100.csv')
eukaryotes_diso_after=pd.read_csv('viruses_long_disorder_after_motif_data_100.csv')

merged_eukaryotes_diso_data=merge_score_data(eukaryotes_diso_before,eukaryotes_diso_after)

print(merged_eukaryotes_diso_data)
merged_eukaryotes_diso_data.to_csv('merged_viruses_long_diso_data_trial.csv', encoding='utf-8', index=False, header= True)