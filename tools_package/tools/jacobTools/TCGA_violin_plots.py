# Imported modules
from importlib.metadata import distribution
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

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

def make_master_df(cancer_list=['BLCA','BRCA', 'COAD', 'GBM', 'KIRC', 'LGG', 'LUAD', 'LUSC', 'OV', 'SARC', 'SKCM', 'UCEC']):
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

    global df_master
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
    if 'LGG' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LGG'], value='Glioma (LGG)')
    if 'LUAD' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LUAD'], value='Lung (LUAD)')
    if 'LUSC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-LUSC'], value='Lung (LUSC)')
    if 'OV' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-OV'], value='Ovarian (OV)')
    if 'SARC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-SARC'], value='Ovarian (SARC)')
    if 'SKCM' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-SKCM'], value='Melanoma (SKCM)')
    if 'UCEC' in cancer_list:
        df_master["Cancer Type"] = df_master["Cancer Type"].replace(['TCGA-UCEC'], value='Uterine (UCEC)')


def make_sub_master_df(aliquot_conc_threshold=0.5):
    #Assign aliquot concentration threshold for all sub master dataframes
    df_master_aliquot = df_master[df_master["Concentration"] >= aliquot_conc_threshold]

    #Make all sub master dataframes
    global df_master_primary
    df_master_primary = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor'])]
    global df_master_blood
    df_master_blood = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Blood Derived Normal'])]
    global df_master_combo
    df_master_combo = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor', 'Blood Derived Normal'])]
    global df_master_all_tumors
    df_master_all_tumors = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Primary Tumor', 'Recurrent Tumor', 'Metastatic'])]
    global df_master_all_normals
    df_master_all_normals = df_master_aliquot.loc[df_master_aliquot['Sample Type'].isin(['Blood Derived Normal', 'Buccal Cell Normal', 'Solid Tissue Normal'])]

    print('Sub master dataframes created!')

def violin_swarm_plot(dataframe, plot_scale=(48,27)):
    sns.set(style="darkgrid")
    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.violinplot(x='Cancer Type', y='Tissue-Specific MMBIR Signatures', saturation=0.75, gridsize=1000, scale="width", data=dataframe, hue_order=["Primary Tumor", "Blood Derived Normal"], hue="Sample Type", inner=None, palette=["cornflowerblue", "tab:red"])  #palette='tab10'
    #Final figure y-axis set below
    #g1.set(ylim=(0, 2000))
    g2 = sns.swarmplot(x='Cancer Type', y='Tissue-Specific MMBIR Signatures', data=dataframe, linewidth = 0.65, alpha=0.5, dodge=True, hue="Sample Type", hue_order=["Primary Tumor", "Blood Derived Normal"], edgecolor="k", s=3, color='k', ax=g1)
    #Final figure size set below
    fig.set_size_inches(16,9)
    plt.legend().remove()
    plt.xlabel("Cancer Type", fontweight="bold", fontsize = 25)   #fontsize = 0
    plt.ylabel("Tissue-Specific MMBIR Signatures", fontweight="bold", fontsize = 25)
    plt.xticks(fontweight="bold", fontsize=12)   #fontsize = 0
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
    fig.set_size_inches(16,9)

    plt.legend().remove()
    plt.xlabel("Cancer Type", fontweight="bold", fontsize = 25)   #fontsize = 0
    plt.ylabel("Tissue-Specific MMBIR Signatures", fontweight="bold", fontsize = 25)
    plt.xticks(fontweight="bold", fontsize=12)   #fontsize = 0
    plt.yticks(fontweight="bold", fontsize=18)
    plt.show()
##################################
def main():
    #List of cancers wanted in master dataframe
    #my_cancer_list=['BLCA','BRCA', 'COAD', 'GBM', 'KIRC', 'LGG', 'LUAD', 'LUSC', 'OV', 'SARC', 'SKCM', 'UCEC']
    make_master_df()
    make_sub_master_df()

if __name__ == '__main__':
    main()

##################################

#List below is just for reference - not used in a function
all_master_dfs = [df_master_primary, df_master_blood, df_master_combo, df_master_all_tumors, df_master_all_normals]

violin_swarm_plot(df_master_primary)
swarm_plot(df_master_primary)



#Modified functions of mmbir_ctrl_vs_tumor_v2.py - might be useful still
def jacob_plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration"):
    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_consolidated = df_consolidated[df_consolidated["cases.0.project.project_id"] == f'TCGA-{cancer}']

    print("Mean and standard deviation of the Filtered_Count:")
    print(df_consolidated.groupby("Sample_Type").mean()["Filtered_Count"])
    print(df_consolidated.groupby("Sample_Type").std()["Filtered_Count"])

    print("Cases with concentration of 0.5:")
    print(df_consolidated[df_consolidated["Concentration"]==.5])

    if hue == "Concentration_bin":
        #bin the concentration into 2 groups for plotting 2 discrete colors based on the concentration
        df_consolidated["Concentration_bin"] = "none"
        df_consolidated["Concentration_bin"] = df_consolidated["Concentration"].apply(lambda x: "low" if x < .5 else "high")

    # set the plot context
    sns.set_context("poster")
    sns.scatterplot(x="Filtered_Count", y="Raw_Count", data=df_consolidated, alpha=0.5, hue=hue, palette="flare")

    # set both x- and y-axis to log scale
    plt.xscale("log")
    plt.yscale("log")

    #rename the x-axis and y-axis
    plt.xlabel("Filtered MMBIR Signature Count")
    plt.ylabel("Raw MMBIR Signature Count")
    plt.title(f"{cancer} Cancer", fontweight="bold")

    plt.show()

def jacob_plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=0.5):


    df_figure = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_figure = df_figure[df_figure["Concentration"] >= min_concentration]

    df_figure = df_figure[df_figure["cases.0.project.project_id"] == f'TCGA-{cancer}']
    
    length=len(df_figure)
    bin_number=int(length/1.5)

    sns.set_context("poster")
    #set color for each sample type
    color_dict = {"Primary Tumor": "tab:blue", "Blood Derived Normal": "tab:orange"}
    sns.histplot(data=df_figure, x="Raw_Count", bins=bin_number, multiple="stack", hue="Sample_Type", palette=color_dict, legend=False) #multiple="stack", hue="Sample_Type",
    
    #make labels bold
    plt.xlabel("Tumor-Specific MMBIR Signatures", fontweight="bold")
    plt.ylabel("Number of Tumor Samples", fontweight="bold")
    
    #add title
    #show the length of the df_figure on the title
    
    plt.title(f"{cancer} Cancer (N={length})", fontweight="bold")

    #don't show the legend
    #plt.legend().remove()

    #set x-axis to go to 1080
    plt.xlim(0,1200)
    plt.ylim(0,12)
    
    #make tick labels bold
    plt.xticks(fontweight="bold", fontsize=20)
    plt.yticks(fontweight="bold", fontsize=20)
    plt.show()
