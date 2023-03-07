# Imported modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Required function from mmbir_ctrl_vs_tumor_v2.py
def annotate_consolidated_results(df_consolidated, df_metadata):

    #join df_metadata and df_consolidated on the file_name column for df_metadata and for Sample_Name for df_consolidated
    df_consolidated=pd.merge(df_consolidated, df_metadata, left_on="Sample_Name", right_on="file_name")

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.diagnoses.0.age_at_diagnosis"] + df_consolidated["cases.0.samples.0.days_to_collection"]
    #divide age_at_collection column by 365 to get the age in years
    df_consolidated["age_at_collection"] = df_consolidated["age_at_collection"]/365

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.demographic.days_to_birth"]
    #divide age_at_collection column by 365 to get the age in years
    df_consolidated["age_at_collection"] = df_consolidated["age_at_collection"]/365 * -1

    #rename the cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration to concentration
    df_consolidated.rename(columns={"cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration": "Concentration"}, inplace=True)

    return df_consolidated

###############
#Function to make all of the dataframes

def make_master_df(cancer_list=['BLCA','BRCA', 'COAD', 'GBM', 'KIRC', 'LGG', 'LUAD', 'LUSC', 'OV', 'SARC', 'SKCM', 'UCEC']): #'LAML'
#Make a master dataframe by concating dfs in cancer_list
#Need to have a consolidated_results anad metadata file for each cancer in the cancer_list
    intial_counter=1
    for cancer in cancer_list:
        if intial_counter == 1:
            print('Starting')
            df_consolidated=pd.read_csv(f"consolidated_results_{cancer}.tsv", sep="\t")
            df_metadata=pd.read_csv(f"TCGA-{cancer}-WXS-BAM-metadata.tsv", sep="\t")
            intial_counter -= 1
            print(df_consolidated)
            print(df_metadata)
        else:
            print(f'Appending {cancer}')
            df_append_consol=pd.read_csv(f"consolidated_results_{cancer}.tsv", sep="\t")
            df_append_meta=pd.read_csv(f"TCGA-{cancer}-WXS-BAM-metadata.tsv", sep="\t")
            print(df_append_consol)
            print(df_append_meta)
            df_consolidated=pd.concat([df_consolidated,df_append_consol])
            df_metadata=pd.concat([df_metadata,df_append_meta])
            print(df_consolidated)
            print(df_metadata)
            print(f'Finshed Appending {cancer}')

    df_master = annotate_consolidated_results(df_consolidated, df_metadata)

    print(f'Master dataframe created!')    

    #Optional column name changes but downstream functions use these new names so be aware
    df_master = df_master.rename(columns={'cases.0.project.project_id':'Cancer Type', 'Raw_Count':'Total MMBIR Signatures', 'Filtered_Count':'Tissue-Specific MMBIR Signatures', 'Sample_Type':'Sample Type'})
    
    #Update cancer type names (not all cancer types included)
    if 'BLCA' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-BLCA'], value='Bladder (BLCA)')
    if 'BRCA' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-BRCA'], value='Breast (BRCA)')
    if 'COAD' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-COAD'], value='Colon (COAD)')
    if 'GBM' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-GBM'], value='Glioblastoma (GBM)')
    if 'KIRC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-KIRC'], value='Kidney (KIRC)')
    if 'LAML' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LAML'], value='Leukemia (LAML)')
    if 'LGG' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LGG'], value='Glioma (LGG)')
    if 'LUAD' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LUAD'], value='Lung (LUAD)')
    if 'LUSC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LUSC'], value='Lung (LUSC)')
    if 'OV' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-OV'], value='Ovarian (OV)')
    if 'PRAD' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-PRAD'], value='Prostate (PRAD)')
    if 'SARC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-SARC'], value='Sarcoma (SARC)')
    if 'SKCM' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-SKCM'], value='Melanoma (SKCM)')
    if 'UCEC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-UCEC'], value='Uterine (UCEC)')
    
    return df_master


def make_sub_master_df(dataframe, aliquot_conc_threshold=0):
    #Assign aliquot concentration threshold for all sub master dataframes
    if aliquot_conc_threshold > 0:
        df_master_aliquot = dataframe[dataframe["Concentration"] >= aliquot_conc_threshold]
    else:
        df_master_aliquot = dataframe

    #Make all sub master dataframes
    df_master_primary = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor'])] #'Additional - New Primary'

    df_master_blood = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Blood Derived Normal'])] #'Primary Blood Derived Cancer - Peripheral Blood'

    df_master_combo = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor', 'Blood Derived Normal'])] #'Primary Blood Derived Cancer - Peripheral Blood'

    df_master_all_tumors = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor', 'Recurrent Tumor', 'Metastatic', 'Additional - New Primary'])]

    df_master_all_normals = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Blood Derived Normal', 'Buccal Cell Normal', 'Solid Tissue Normal'])]

    print('Sub master dataframes created!')

    return df_master_aliquot, df_master_primary, df_master_blood, df_master_combo, df_master_all_tumors, df_master_all_normals

def violin_swarm_plot(dataframe, plot_scale=(48,27)):
    sns.set(style="darkgrid")
    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.violinplot(x='Cancer Type', y='Tissue-Specific MMBIR Signatures', saturation=0.75, gridsize=1000, scale="width", data=dataframe, hue_order=["Primary Tumor",'Blood Derived Normal'], hue="Sample Type", inner=None, palette=["#E5EEFB", "mistyrose"])  #palette='tab10', 'Blood Derived Normal'
    #Final figure y-axis set below
    g1.set(ylim=(0, 7000))
    
    g2 = sns.swarmplot(x='Cancer Type', y='Tissue-Specific MMBIR Signatures', data=dataframe, linewidth = 0.5, alpha=0.5, dodge=True, hue="Sample Type", hue_order=["Primary Tumor", 'Blood Derived Normal'], edgecolor="grey", s=3, palette=["royalblue", "tab:red"], ax=g1) #, 'Blood Derived Normal'
    #Final figure size set below
    #fig.set_size_inches(16,9)
    plt.legend().remove()
    #plt.xlabel("Cancer Type", fontweight="bold", fontsize = 25)   #fontsize = 0
    plt.ylabel("Tissue-Specific MMBIR Signatures", fontweight="bold", fontsize = 25)
    plt.xticks(fontweight="bold", fontsize=20)   #fontsize = 0
    plt.yticks(fontweight="bold", fontsize=18)
    plt.show()

def swarm_plot(dataframe, plot_scale=(48,27)):
    sns.set(style="darkgrid")
    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.swarmplot(x='Cancer Type', y='Tissue-Specific MMBIR Signatures', data=dataframe, linewidth = 0.65, alpha=0.5, dodge=True, hue="Sample Type", hue_order=["Primary Tumor", "Blood Derived Normal"], edgecolor="k", s=3, palette=["cornflowerblue", "tab:red"])  #palette='tab10'
    
    #Final figure y-axis set below
    #g1.set(ylim=(0, 2000))

    #Final figure size set below
    #fig.set_size_inches(16,9)

    plt.legend().remove()
    plt.xlabel("Cancer Type", fontweight="bold", fontsize = 25)   #fontsize = 0
    plt.ylabel("Tissue-Specific MMBIR Signatures", fontweight="bold", fontsize = 25)
    plt.xticks(fontweight="bold", fontsize=12)   #fontsize = 0
    plt.yticks(fontweight="bold", fontsize=18)
    plt.show()

def calculate_n(df):
    analyzed_cancers = set(df["Cancer Type"].values.tolist())
    sample_types = set((df["Sample Type"].values.tolist())[::-1])
    for cancer in analyzed_cancers:
        print(f"For {cancer} cancer:")
        df_temp = df.loc[df["Cancer Type"]==cancer]
        for location in sample_types:
            n_value = df_temp["Sample Type"].value_counts()[location]
            print(f"\t{location} n value = {n_value}")

def calculate_mean(df):
    analyzed_cancers = set(df["Cancer Type"].values.tolist())
    sample_types = set((df["Sample Type"].values.tolist())[::-1])
    for cancer in analyzed_cancers:
        print(f"For {cancer} cancer:")
        df_temp = df.loc[df["Cancer Type"]==cancer]
        for location in sample_types:
            df_temp2 = df_temp.loc[df_temp["Sample Type"].isin([location])]
            mean = round(df_temp2["Tissue-Specific MMBIR Signatures"].mean())
            print(f"\t{location} average MMBIR Signatures = {mean}")

##################################
def main():
    #List of cancers wanted in master dataframe
    #my_cancer_list=['BRCA', 'COAD', 'GBM', 'KIRC', 'LUAD', 'OV', 'UCEC'] #'BLCA', 'LAML', 'LGG',  'LUSC',  'SARC', 'SKCM',
    df_master = make_master_df()
    df_tuple = make_sub_master_df(df_master, 0.5) #Value in fucntion is the aliquot thershold cutoff
    df_master_aliquot = df_tuple[0]
    df_master_primary = df_tuple[1]
    df_master_blood = df_tuple[2]
    df_master_combo = df_tuple[3]
    df_master_all_tumors = df_tuple[4]
    df_master_all_normals = df_tuple[5]

    #all_master_dfs = [df_master_primary, df_master_blood, df_master_combo, df_master_all_tumors, df_master_all_normals]

    swarm_plot(df_master_combo)
    violin_swarm_plot(df_master_combo)
    calculate_mean(df_master_combo)

if __name__ == '__main__':
    main()

##################################