import pandas as pd
import logging
import scipy
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
from numpy.lib import scimath


expression_df = pd.read_csv("expression_data_with_CaseID.tsv", sep="\t")
consolidated_results_path = "consolidated_results.tsv"
MMBIR_THRESHOLD=300


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


threshold_mmb_cases_df = getCasesAboveMMBThreshold(consolidated_results_path, MMBIR_THRESHOLD)

# get the IDs of the cases that are high mmbir
high_mmbir_cases=threshold_mmb_cases_df["Case_ID_"].values.tolist()

# Mark the cases that are above the MMBIR threshold as high in the expression dataframe
expression_df["high_mmbir"]=expression_df["case_id"].apply(lambda x: 1 if x in high_mmbir_cases else 0)


# get the list of all columns (transcripts) in the expression dataframe
expression_transcripts = expression_df.columns.values.tolist()[:-3]

#expression_transcripts=["ROBO3"]

expression_dict = {}
for transcript in expression_transcripts:

    # focus on the expression dataframe for the transcript
    transcript_df = expression_df[[transcript, "high_mmbir"]].copy()

    # take the log2 of the transcript expression values and add 1 to avoid log2(0)
    transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+1))

    # create violin plots for the transcript_df in seaborn
    '''sns.violinplot(data=transcript_df, x="high_mmbir", y=transcript)
    plt.show()
    sns.swarmplot(data=transcript_df, x="high_mmbir", y=transcript)
    plt.show()'''

    # split the expression dataframe into high and low mmbir
    high_mmbir_transcript_df = transcript_df[transcript_df["high_mmbir"] == 1]
    low_mmbir_transcript_df = transcript_df[transcript_df["high_mmbir"] == 0]

    # get the mean expression values for the high and low mmbir cases
    high_mmbir_mean = high_mmbir_transcript_df[f"{transcript}_log2"].mean()
    low_mmbir_mean = low_mmbir_transcript_df[f"{transcript}_log2"].mean()

    high_mmbir_transcript_df = high_mmbir_transcript_df[f"{transcript}_log2"]
    #print(f"High: {high_mmbir_transcript_df.describe()}")
    low_mmbir_transcript_df = low_mmbir_transcript_df[f"{transcript}_log2"]
    #print(f"Low: {low_mmbir_transcript_df.describe()}")


    # Perform the Mann-Whitney U rank test on the expression dataframes
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
    try:
        fold_change = high_mmbir_mean/low_mmbir_mean
    except ZeroDivisionError:
        fold_change = "NA"

    result = scipy.stats.mannwhitneyu(high_mmbir_transcript_df, low_mmbir_transcript_df, nan_policy='omit')
    print(f"{transcript}, U: {result[0]}, p-value: {result[1]}, fold-change: {fold_change}")

    # store the results and fold change in a dictionary
    expression_dict[transcript] = {"U": result[0], "p-value": result[1], "fold-change": fold_change}



# sort the dictionary by the p-value of the Mann-Whitney U test and output the results to a list
sorted_expression_list = sorted(expression_dict.items(), key=lambda x: x[1]["p-value"])

#save the results to a file
with open("mannwhitneyu_results_new_300.tsv", "w") as f:
    f.write("Transcript\tU\tp-value\tfold-change\n")
    for transcript, result in sorted_expression_list:
        f.write(f"{transcript}\t{result['U']}\t{result['p-value']}\t{result['fold-change']}\n")
