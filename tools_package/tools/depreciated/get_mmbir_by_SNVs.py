# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import numpy
import os
import logging

snvs_metadata="/Users/twarowski/MMBIR_Databases/TCGA/TCGA-BRCA-WXS-MAF-metadata.tsv"
snvs_loc="/nfsscratch/twarowski/TCGA-BRCA/Masked_MAF"
df_metadata = pd.read_csv(snvs_metadata, sep="\t")

output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
print(f"cases exited with a code {output}")
cases_file = open("cases.txt", 'r')
cases = cases_file.readlines()
cases_file.close()

consolidated_results_path = "consolidated_results.tsv"
MMBIR_THRESHOLD=200


def extractAffectedGeneNames(maf_path):

    try:
        df_masked_snvs=pd.read_csv(maf_path, sep="\t", comment='#')
    except OSError:
        logging.info(f"Couldnt open the MAF file: {maf_path} Is it there?")
        return None
    try:
        # remove silent mutations
        df_masked_snvs = df_masked_snvs[df_masked_snvs["Variant_Classification"] != "Silent"]
        return_gene_list = df_masked_snvs["Hugo_Symbol"].values.tolist()
        return return_gene_list
    except ValueError:
        logging.info(f"Couldn't parse the MAF file: {maf_path}")
        return None


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
###################################################################################



threshold_mmb_cases_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD)

high_mmbir_cases=0
high_mmbir_snv_genes=[]
low_mmbir_cases=0
low_mmbir_snv_genes=[]

for caseID in cases:

    caseID=caseID.strip()
    case_snvs = df_metadata[df_metadata["cases.0.case_id"] == caseID]
    case_files = case_snvs[["file_name", "id"]].values.tolist()

    #get all the snv affected genes (for now just pick the first maf file comming in)
    for file_name, id in case_files:
        maf_path = f"{caseID}/{file_name}"
        maf_path = maf_path.strip(".gz")
        maf_genes = extractAffectedGeneNames(maf_path)


    #check if case is high or low mmbir
    if caseID in threshold_mmb_cases_df["Case_ID_"].values.tolist():

        #if yes -> HIGH MMBIR, count as high mmbir sample, add to high mmbir cohort
        high_mmbir_cases+=1
        high_mmbir_snv_genes.extend(maf_genes)
    else:

        #if no -> LOW MMBIR, count as low mmbir sample, add to low mmbir cohort
        low_mmbir_cases+=1
        low_mmbir_snv_genes.extend(maf_genes)



fdist_HIGH=dict(zip(*numpy.unique(high_mmbir_snv_genes, return_counts=True)))
fdist_HIGH = {k: v / high_mmbir_cases for k, v in fdist_HIGH.items()}
sorted_fdist_HIGH = sorted(fdist_HIGH.items(), key=lambda x: x[1], reverse=True)

fdist_LOW=dict(zip(*numpy.unique(low_mmbir_snv_genes, return_counts=True)))
fdist_LOW = {k: v / low_mmbir_cases for k, v in fdist_LOW.items()}
sorted_fdist_LOW = sorted(fdist_LOW.items(), key=lambda x: x[1], reverse=True)

fdist_all=dict(zip(*numpy.unique(high_mmbir_snv_genes+low_mmbir_snv_genes, return_counts=True)))
fdist_all = {k: v / (high_mmbir_cases+low_mmbir_cases) for k, v in fdist_all.items()}
sorted_fdist_all = sorted(fdist_all.items(), key=lambda x: x[1], reverse=True)


print(f"found HIGH: {high_mmbir_cases} and LOW: {low_mmbir_cases} cases")
print("HIGH")
print("Hugo,high-set,low-set,all-set")
for i in sorted_fdist_HIGH:
    try:
        compl_low=fdist_LOW[i[0]]
        compl_all=fdist_all[i[0]]
    except:
        logging.info(f"No corresponding record found for {i[0]} in the LOW-MMBIR Set")
        compl_low="NULL"
    print(f"{i[0]},{i[1]},{compl_low},{compl_all}")


print("LOW")
for i in sorted_fdist_LOW:
    try:
        compl_high=fdist_HIGH[i[0]]
        compl_all=fdist_all[i[0]]
    except:
        logging.info(f"No corresponding record found for {i[0]} in the HIGH-MMBIR Set")
        compl_high="NULL"
    print(f"{i[0]},{compl_high},{i[1]},{compl_all}")
