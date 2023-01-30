import pandas as pd
import logging
import scipy
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
from numpy.lib import scimath
import cancer_config as cfg
from tools import getCasesAboveMMBThreshold

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

consolidated_results_path = cfg.settings["consolidated_results_path"]
expression_df_path = cfg.settings["expression_df_path"]

expression_df = pd.read_csv(expression_df_path, sep="\t")
df_sample_metadata = pd.read_csv(metadata_location, sep="\t")


#threshold_mmb_cases_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD)

threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_HIGH, min_concentration=0)
threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD_LOW, below=True, min_concentration=0)


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

expression_transcripts=["MSH5"]

expression_dict = {}
for transcript in expression_transcripts:

    # focus on the expression dataframe for the transcript
    transcript_df = expression_df[[transcript, "high_mmbir"]].copy()

    # take the log2 of the transcript expression values and add 0.1 to avoid log2(0)
    transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+0.1))

    # create violin plots for the transcript_df in seaborn
    sns.violinplot(data=transcript_df, x="high_mmbir", y=transcript)
    plt.show()
    sns.swarmplot(data=transcript_df, x="high_mmbir", y=transcript)
    plt.show()

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


    #print(f"{transcript}, t: {t_stat}, p-value: {p_value}, fold-change: {fold_change}")

    # store the results and fold change in a dictionary
    expression_dict[transcript] = {"t": t_stat, "p-value": p_value, "fold-change": fold_change}



# sort the dictionary by the two-sample t-test and output the results to a list
sorted_expression_list = sorted(expression_dict.items(), key=lambda x: x[1]["p-value"])

#save the results to a file
with open("t_test_200_300_result_highconc.tsv", "w") as f:
    f.write("Transcript\tt\tp-value\tfold-change\n")
    for transcript, result in sorted_expression_list:
        f.write(f"{transcript}\t{result['t']}\t{result['p-value']}\t{result['fold-change']}\n")

