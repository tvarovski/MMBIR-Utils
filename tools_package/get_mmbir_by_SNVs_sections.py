import pandas as pd
import numpy
import os
import logging
import cancer_config as cfg
from tools import extractAffectedGeneNames, getCasesAboveMMBThreshold

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

consolidated_results_path = cfg.settings["consolidated_results_path"]

snvs_metadata=f"/Users/{username}/MMBIR_Databases/TCGA/TCGA-{cancer}-WXS-MAF-metadata.tsv"
snvs_loc=f"/nfsscratch/{username}/TCGA-{cancer}/Masked_MAF"

df_metadata = pd.read_csv(snvs_metadata, sep="\t")
df_sample_metadata = pd.read_csv(metadata_location, sep="\t")


output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
print(f"cases exited with a code {output}")
cases_file = open("cases.txt", 'r')
cases = cases_file.readlines()
cases_file.close()

###################################################################################

threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_HIGH, min_concentration=0)
threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_LOW, below=True, , min_concentration=0)

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
    if caseID in threshold_mmb_cases_high_df["Case_ID_"].values.tolist():

        #if yes -> HIGH MMBIR, count as high mmbir sample, add to high mmbir cohort
        high_mmbir_cases+=1
        high_mmbir_snv_genes.extend(maf_genes)
    elif caseID in threshold_mmb_cases_low_df["Case_ID_"].values.tolist():

        #if no -> LOW MMBIR, count as low mmbir sample, add to low mmbir cohort
        low_mmbir_cases+=1
        low_mmbir_snv_genes.extend(maf_genes)
    else:
        #if neither -> set is ambiguous, do nothing
        pass



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

