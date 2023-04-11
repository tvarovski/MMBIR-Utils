# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import logging
import scipy
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
from numpy.lib import scimath
import statsmodels.stats.multitest as smm
import cancer_config as cfg
from tools import getCasesAboveMMBThreshold

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

expression_df_path = cfg.settings["expression_df_path"]
consolidated_results_path = cfg.settings["consolidated_results_path"]

df_sample_metadata = pd.read_csv(metadata_location, sep="\t")
expression_df = pd.read_csv(expression_df_path, sep="\t")


#threshold_mmb_cases_high_df = df_consolidated=pd.read_csv("case_ids_rad52_SNP.txt", sep="\t")
#threshold_mmb_cases_low_df = df_consolidated=pd.read_csv("case_ids_rad52_wtSNP.txt", sep="\t")

#threshold_mmb_cases_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD)

threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_HIGH, min_concentration=0)
threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_LOW, below=True, min_concentration=0)


# get the IDs of the cases that are high mmbir
high_mmbir_cases=threshold_mmb_cases_high_df["Case_ID_"].values.tolist()
low_mmbir_cases=threshold_mmb_cases_low_df["Case_ID_"].values.tolist()


# Mark the cases that are above the MMBIR threshold as high in the expression dataframe
expression_df["high_mmbir"]="none"
expression_df.loc[expression_df["case_id"].isin(high_mmbir_cases), "high_mmbir"]="high"
expression_df.loc[expression_df["case_id"].isin(low_mmbir_cases), "high_mmbir"]="low"

len_high=len(expression_df[expression_df["high_mmbir"] == "high"])
len_low=len(expression_df[expression_df["high_mmbir"] == "low"])
print(f"high cases: {len_high}, low cases: {len_low}")


# get the list of all columns (transcripts) in the expression dataframe
expression_transcripts = expression_df.columns.values.tolist()[:-3]

#expression_transcripts=["RAD52", "BRCA1", "BRCA2", "TP53"]

expression_dict = {}
for transcript in expression_transcripts:

    # focus on the expression dataframe for the transcript
    transcript_df = expression_df[[transcript, "high_mmbir"]].copy()

    # take the log2 of the transcript expression values and add 0.1 to avoid log2(0)
    transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+0.1))

    # create violin plots for the transcript_df in seaborn
    #sns.violinplot(data=transcript_df, x="high_mmbir", y=transcript)
    #plt.show()
    #sns.swarmplot(data=transcript_df, x="high_mmbir", y=transcript)
    #plt.show()

    # split the expression dataframe into high and low mmbir
    high_mmbir_transcript_df = transcript_df[transcript_df["high_mmbir"] == "high"]
    low_mmbir_transcript_df = transcript_df[transcript_df["high_mmbir"] == "low"]

    # get the mean expression values for the high and low mmbir cases
    high_mmbir_mean = high_mmbir_transcript_df[f"{transcript}_log2"].mean()
    low_mmbir_mean = low_mmbir_transcript_df[f"{transcript}_log2"].mean()

    high_mmbir_transcript_df = high_mmbir_transcript_df[f"{transcript}_log2"]
    #print(f"High: {high_mmbir_transcript_df.describe()}")
    low_mmbir_transcript_df = low_mmbir_transcript_df[f"{transcript}_log2"]
    #print(f"Low: {low_mmbir_transcript_df.describe()}")


    # Perform t-test to see if the mean expression values for the high and low mmbir cases are different
    t_stat, p_value = scipy.stats.ttest_ind(high_mmbir_transcript_df, low_mmbir_transcript_df, nan_policy='omit')
    try:
        fold_change = high_mmbir_mean/low_mmbir_mean
    except ZeroDivisionError:
        fold_change = "NA"


    print(f"{transcript}, t: {t_stat}, p-value: {p_value}, fold-change: {fold_change}")

    # store the results and fold change in a dictionary
    expression_dict[transcript] = {"t": t_stat, "p-value": p_value, "fold-change": fold_change, "high_mmbir_mean": high_mmbir_mean, "low_mmbir_mean": low_mmbir_mean}


# sort the dictionary by the two-sample t-test and output the results to a list
sorted_expression_list = sorted(expression_dict.items(), key=lambda x: x[1]["p-value"])

OUTPUT_NAME="t_test_High250_Low150.tsv"

#save the results to a file
with open(OUTPUT_NAME, "w") as f:
    f.write("Transcript\tt\tp-value\tfold-change\tmmbir_high_mean\tmmbir_low_mean\n")
    for transcript, result in sorted_expression_list:
        f.write(f"{transcript}\t{result['t']}\t{result['p-value']}\t{result['fold-change']}\t{result['high_mmbir_mean']}\t{result['low_mmbir_mean']}\n")

expression_results_df=pd.read_csv(OUTPUT_NAME, sep="\t")
expression_results_df=expression_results_df.sort_values(by=['p-value'], ascending=True)
expression_results_df.reset_index(drop=True, inplace=True)
expression_results_df['p-value_corrected'] = smm.multipletests(expression_results_df['p-value'], method='fdr_bh')[1]

expression_results_df.to_csv(OUTPUT_NAME, sep="\t", index=False)
