import pandas as pd
import os
import logging
import cancer_config as cfg
from tools import getCasesAboveMMBThreshold, returnGenes, countGenes

#returnGenes needs update in tools_base

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

raw_dir = cfg.setting["outputs_path_raw"]
filtered_dir = cfg.setting["outputs_path_filtered"]

output_file_name="../gene_frequencies_mmbir_sections_highconc_test.csv"

#os.chdir("/Shared/malkova_lab/Jacob/TCGA_Glioblastoma_Project/")

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
threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_LOW, min_concentration, below=True)

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
    genes = returnGenes(file, high_mmbir, exones_only, ref_homology_check_filter)
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

# sort the dictionaries by frequency of occurence of genes in the raw and filtered files
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
df_frequency.to_csv(output_file_name)