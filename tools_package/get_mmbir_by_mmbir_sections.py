import pandas as pd
import os
import logging
import cancer_config as cfg

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

consolidated_results_path = cfg.settings["consolidated_results_path"]

df_metadata = pd.read_csv(metadata_location, sep="\t")

ref_complexity_filter = cfg.settings["ref_complexity_filter"]
bir_complexity_filter =  cfg.settings["bir_complexity_filter"]
ref_homology_check_filter = cfg.settings["ref_homology_check_filter"]
bir_homology_check_filter = cfg.settings["bir_homology_check_filter"]
exones_only = cfg.settings["exones_only"]

min_concentration = 0.0

#### Needs Refactoring with Up ####
exones_only = True
homology_ref = False
homology_bir = False
complexity_ref = True
complexity_bir = True

raw_dir="outputs/raw"
filtered_dir="outputs/filtered"

def getCasesAboveMMBThreshold(consolidated_results_path, min_MMBIR_events, below=False, min_concentration=0):

    df_consolidated=pd.read_csv(consolidated_results_path, sep="\t")

    df_consolidated=pd.merge(df_consolidated, df_sample_metadata, left_on="Sample_Name", right_on="file_name")

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.diagnoses.0.age_at_diagnosis"] + df_consolidated["cases.0.samples.0.days_to_collection"]
    df_consolidated.rename(columns={"cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration": "Concentration"}, inplace=True)
    df_consolidated=df_consolidated[df_consolidated["Concentration"] >= min_concentration]


    agg_dict={"Raw_Count": ['min', 'max'],
              "Filtered_Count": ['min', 'max']}
    df_agg = df_consolidated.groupby("Case_ID").agg(agg_dict).reset_index()

    df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]

    if below:
        df_agg = df_agg[df_agg["Filtered_Count_max"] < min_MMBIR_events]
    else:
        df_agg = df_agg[df_agg["Filtered_Count_max"] >= min_MMBIR_events]


    logging.info(f"The length is: {len(df_agg)}")

    df_agg = df_agg.sort_values("Filtered_Count_max",ascending=False)

    return df_agg

def returnGenes(file, high_mmbir, exones_only=False, homology_ref=False, homology_bir=False, complexity_ref=False, complexity_bir=False):

    df = pd.read_csv(file, sep='\t')

    if exones_only:
        df = df[df["exones"].str.len() != 2]
        #print(df)
    if homology_ref:
        df = df[df["homology_check_ref"] == True]
    if homology_bir:
        df = df[df["homology_check_bir"] == True]
    if complexity_ref:
        df = df[df["ref_complexity_fail"] == False]
    if complexity_bir:
        df = df[df["bir_complexity_fail"] == False]

    # read genes column of the df and return list of genes
    genes = df["genes"].tolist()
    # get unique genes
    genes = list(set(genes))

    #if mmbir is high, return each gene with a high flag
    if high_mmbir == True:
        for index, gene in enumerate(genes):
            genes[index] = (gene, True)
    
    #if mmbir is low, return each gene with a low flag
    elif high_mmbir == False:
        for index, gene in enumerate(genes):
            genes[index] = (gene, False)
    
    else:
        raise ValueError("high_mmbir must be True or False")

    return genes

def countGenes(genes):
    gene_count_high = {}
    gene_count_low = {}
    for gene, high_mmbir in genes:
        if high_mmbir == True:
            if gene in gene_count_high:
                gene_count_high[gene] += 1
            else:
                gene_count_high[gene] = 1

        if high_mmbir == False:
            if gene in gene_count_low:
                gene_count_low[gene] += 1
            else:
                gene_count_low[gene] = 1
    return gene_count_high, gene_count_low


os.chdir("/Shared/malkova_lab/Jacob/TCGA_Glioblastoma_Project/")

#get the current working directory
cwd = os.getcwd()

# add the raw_dir to the current working directory
os.chdir(cwd+"/"+raw_dir)

#list all files in the raw_dir ending with .txt
raw_files = [f for f in os.listdir() if f.endswith(".txt")]

#displaythe number of files in the raw_dir
print("There are {} files in the raw_dir".format(len(raw_files)))

os.chdir(cwd)

os.chdir(cwd+"/"+filtered_dir)
filtered_files = [f for f in os.listdir() if f.endswith(".txt")]
os.chdir(cwd)


# get the IDs of the cases that are high mmbir
threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_HIGH, min_concentration)
threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_LOW, below=True, min_concentration)



high_mmbir_cases = threshold_mmb_cases_high_df["Case_ID_"].values.tolist()
low_mmbir_cases = threshold_mmb_cases_low_df["Case_ID_"].values.tolist()

all_high_mmbir_files = df_metadata[df_metadata["cases.0.case_id"].isin(high_mmbir_cases)].file_name.values.tolist()
all_low_mmbir_files = df_metadata[df_metadata["cases.0.case_id"].isin(low_mmbir_cases)].file_name.values.tolist()

# get the number of files associated with high mmbir cases
high_mmbir_files = [f for f in raw_files if f"{os.path.splitext(f)[0][:-4]}.bam" in all_high_mmbir_files]
low_mmbir_files = [f for f in raw_files if f"{os.path.splitext(f)[0][:-4]}.bam" in all_low_mmbir_files]

print(len(high_mmbir_files))

print(f"There are {len(high_mmbir_files)} high mmbir files")
print(f"There are {len(low_mmbir_files)} low mmbir files")

raw_genes=[]
filtered_genes=[]

for file in raw_files:

    # get the file name without the extension
    file_name = os.path.splitext(file)[0]

    file_name = file_name[:-4]+".bam"
    print(file_name)


    # get the case id of the file by referencing metadata file
    case_id = df_metadata[df_metadata["file_name"] == file_name]["cases.0.case_id"].values[0]

    if case_id in high_mmbir_cases:
        high_mmbir=True
    elif case_id in low_mmbir_cases:
        high_mmbir=False
    else:
        high_mmbir="NA"
        continue

    file=raw_dir+"/"+file
    genes = returnGenes(file, high_mmbir, exones_only, homology_ref)
    raw_genes.extend(genes)

for file in filtered_files:

    # get the file name without the extension
    file_name = os.path.splitext(file)[0]
    file_name = file_name[:-25]+".bam"
    print(file_name)
    # get the case id of the file by referencing metadata file
    case_id = df_metadata[df_metadata["file_name"] == file_name]["cases.0.case_id"].values[0]

    if case_id in high_mmbir_cases:
        high_mmbir=True
    elif case_id in low_mmbir_cases:
        high_mmbir=False
    else:
        high_mmbir="NA"
        continue

    file=filtered_dir+"/"+file
    genes = returnGenes(file, high_mmbir, exones_only)
    filtered_genes.extend(genes)

# count the number of occurences of each gene in the raw and filtered files and return a dictionary

gene_freq_raw_high, gene_freq_raw_low = countGenes(raw_genes)
gene_freq_filtered_high, gene_freq_filtered_low = countGenes(filtered_genes)

# sort the dictionaries by frequency of occurence of genes
# in the raw and filtered files
# and return a list of tuples sorted by frequency of occurence
raw_sorted_high = sorted(gene_freq_raw_high.items(), key=lambda x: x[1], reverse=True)
raw_sorted_low = sorted(gene_freq_raw_low.items(), key=lambda x: x[1], reverse=True)
filtered_sorted_high = sorted(gene_freq_filtered_high.items(), key=lambda x: x[1], reverse=True)
filtered_sorted_low = sorted(gene_freq_filtered_low.items(), key=lambda x: x[1], reverse=True)

# divide the frequency of occurence of genes in the raw and filtered files by length of the raw files
# and return a list of tuples sorted by frequency of occurence
raw_sorted_high = [(x[0], x[1]/len(high_mmbir_files)) for x in raw_sorted_high]
raw_sorted_low = [(x[0], x[1]/len(low_mmbir_files)) for x in raw_sorted_low]
filtered_sorted_high = [(x[0], x[1]/len(high_mmbir_files)) for x in filtered_sorted_high]
filtered_sorted_low = [(x[0], x[1]/len(low_mmbir_files)) for x in filtered_sorted_low]


# print the top_N_genes in the raw and filtered files
top_n_genes=10

print(f"Top {top_n_genes} genes in the raw+high data set:")
for i in range(top_n_genes):
    print(raw_sorted_high[i])

print(f"Top {top_n_genes} genes in the raw+low data set:")
for i in range(top_n_genes):
    print(raw_sorted_low[i])

print(f"Top {top_n_genes} genes in the filtered+high data set:")
for i in range(top_n_genes):
    print(filtered_sorted_high[i])

print(f"Top {top_n_genes} genes in the filtered+low data set:")
for i in range(top_n_genes):
    print(filtered_sorted_low[i])


# create a dictionary with gene as key and the frequency of occurence as value
gene_freq_raw_high = dict(raw_sorted_high)
gene_freq_raw_low = dict(raw_sorted_low)
gene_freq_filtered_high = dict(filtered_sorted_high)
gene_freq_filtered_low = dict(filtered_sorted_low)

# get the list of keys in the raw high and low and filtered high and low dictionaries
raw_high_genes = list(gene_freq_raw_high.keys())
raw_low_genes = list(gene_freq_raw_low.keys())
filtered_high_genes = list(gene_freq_filtered_high.keys())
filtered_low_genes = list(gene_freq_filtered_low.keys())

#combine the lists of keys
all_genes = raw_high_genes + raw_low_genes + filtered_high_genes + filtered_low_genes

#get unique list of genes
unique_genes = list(set(all_genes))

# create a nested dictionary with gene as the first key, dataset as the second key
# and the frequency of occurence as value

frequency_dict = {}
#iterate through genes
for gene in unique_genes:
    #retrieve the frequency of occurence of the gene in the raw high and low dictionaries
    raw_high_freq = gene_freq_raw_high.get(gene, 0)
    raw_low_freq = gene_freq_raw_low.get(gene, 0)
    #retrieve the frequency of occurence of the gene in the filtered high and low dictionaries
    filtered_high_freq = gene_freq_filtered_high.get(gene, 0)
    filtered_low_freq = gene_freq_filtered_low.get(gene, 0)

    frequency_dict[gene] = {"raw_high": raw_high_freq,
                             "raw_low": raw_low_freq,
                             "filtered_high": filtered_high_freq,
                             "filtered_low": filtered_low_freq}

# save the frequency dictionary as a dataframe
df_frequency = pd.DataFrame.from_dict(frequency_dict, orient="index")

#order the dataframe by filtered_high frequency
df_frequency = df_frequency.sort_values(by="filtered_high", ascending=False)

# save the dataframe as a csv
df_frequency.to_csv("../gene_frequencies_mmbir_sections_highconc_test.csv")


