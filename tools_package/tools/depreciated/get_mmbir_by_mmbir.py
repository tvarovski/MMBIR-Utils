# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import os
import logging

exone_metadata="/Users/twarowski/MMBIR_Databases/TCGA/TCGA-BRCA-WXS-BAM-metadata.tsv"
# read metadata file
df_metadata = pd.read_csv(exone_metadata, sep="\t")

consolidated_results_path = "consolidated_results.tsv"
MMBIR_THRESHOLD=150

raw_dir="outputs/raw"
filtered_dir="outputs/filtered"

def getCasesAboveMMBThreshold(consolidated_results_path, min_MMBIR_events):

    df_consolidated=pd.read_csv(consolidated_results_path, sep="\t")


    agg_dict={"Raw_Count": ['min', 'max'],
              "Filtered_Count": ['min', 'max']}
    df_agg = df_consolidated.groupby("Case_ID").agg(agg_dict).reset_index()

    df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]

    df_agg = df_agg[df_agg["Filtered_Count_max"] > min_MMBIR_events]

    logging.info(f"The length is: {len(df_agg)}")

    df_agg = df_agg.sort_values("Filtered_Count_max",ascending=False)

    return df_agg

def returnGenes(file, high_mmbir):

    df = pd.read_csv(file, sep='\t')

    # read genes column of the df and return list of genes
    genes = df["genes"].tolist()
    # get unique genes
    genes = list(set(genes))

    #if mmbir is high, return each gene with a high flag
    if high_mmbir:
        for index, gene in enumerate(genes):
            genes[index] = (gene, True)
    else:
        for index, gene in enumerate(genes):
            genes[index] = (gene, False)

    return genes

def countGenes(genes):
    gene_count_high = {}
    gene_count_low = {}
    for gene, high_mmbir in genes:
        if high_mmbir:
            if gene in gene_count_high:
                gene_count_high[gene] += 1
            else:
                gene_count_high[gene] = 1

        if not high_mmbir:
            if gene in gene_count_low:
                gene_count_low[gene] += 1
            else:
                gene_count_low[gene] = 1
    return gene_count_high, gene_count_low


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
threshold_mmb_cases_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD)
high_mmbir_cases = threshold_mmb_cases_df["Case_ID_"].values.tolist()
all_high_mmbir_files = df_metadata[df_metadata["cases.0.case_id"].isin(high_mmbir_cases)].file_name.values.tolist()

# get the number of files associated with high mmbir cases
high_mmbir_files = [f for f in raw_files if f"{os.path.splitext(f)[0][:-4]}.bam" in all_high_mmbir_files]
print(len(high_mmbir_files))

print(f"There are {len(high_mmbir_files)} high mmbir files")
print(f"There are {len(raw_files) - len(high_mmbir_files)} low mmbir files")

low_mmbir_files_count = len(raw_files) - len(high_mmbir_files)


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
    else:
        high_mmbir=False

    file=raw_dir+"/"+file
    genes = returnGenes(file, high_mmbir)
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
    else:
        high_mmbir=False

    file=filtered_dir+"/"+file
    genes = returnGenes(file, high_mmbir)
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
raw_sorted_low = [(x[0], x[1]/low_mmbir_files_count) for x in raw_sorted_low]
filtered_sorted_high = [(x[0], x[1]/len(high_mmbir_files)) for x in filtered_sorted_high]
filtered_sorted_low = [(x[0], x[1]/low_mmbir_files_count) for x in filtered_sorted_low]


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
df_frequency.to_csv("gene_frequencies_mmbir.csv")
