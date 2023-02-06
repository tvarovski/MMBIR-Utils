import pandas as pd
import scipy
import numpy
import statsmodels.stats.multitest
import cancer_config as cfg
from tools import getCasesAboveMMBThreshold


cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

consolidated_results_path = cfg.settings["consolidated_results_path"]
expression_df_path = f"expression_data_{cancer}.pickle"

expression_df = pd.read_pickle(expression_df_path)
df_sample_metadata = pd.read_csv(metadata_location, sep="\t")

output_name = f"ttest_results_{cancer}.pickle"


min_concentration=0
output_name = f"ttest_results_{cancer}_conc{min_concentration}.tsv"

def performTTest(expression_df, df_sample_metadata, min_concentration=0):

    threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_HIGH, min_concentration=min_concentration)
    threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_LOW, below=True, min_concentration=min_concentration)


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

    #expression_transcripts=["MSH5"]

    expression_dict = {}
    for transcript in expression_transcripts:

        # focus on the expression dataframe for the transcript
        transcript_df = expression_df[[transcript, "high_mmbir"]].copy()

        # take the log2 of the transcript expression values; add 0.1 to avoid log2(0)
        try:
            #set transcript_df[transcript] column to float
            transcript_df[transcript] = transcript_df[transcript].astype(float)
            transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+0.1))
        
        except:
            
            #print column names
            print(transcript_df.columns.values.tolist())

            #if there are 2 of the same column name, drop the one whose values are all 0s or NaNs
            if transcript_df.columns.values.tolist().count(transcript) > 1:

                print(f"dropping one of the columns of {transcript}")
                transcript_df = transcript_df.loc[:,~transcript_df.columns.duplicated()]

                transcript_df[transcript] = transcript_df[transcript].astype(float)
                transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+0.1))
            else:
                print(f"couldn't drop a column. Skipping {transcript}")
                continue

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


        #print(f"{transcript}, t: {t_stat}, p-value: {p_value}, fold-change: {fold_change}")

        # store the results and fold change in a dictionary
        expression_dict[transcript] = {"t": t_stat, "p-value": p_value, "fold-change": fold_change}



    # create a dataframe from the dictionary where the index is the transcript name and the columns are the results
    expression_df = pd.DataFrame.from_dict(expression_dict, orient="index")

    # sort the dataframe by the two-sample t-test, lowest p-value first, and output the results to a file
    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(output_name, sep="\t")
    return expression_df


def performBenjaminiHochbergCorrection(expression_df, output_name, *min_p_value):
    # get the list of transcripts that are significantly differentially expressed,
    #perform a benjamini-hochberg correction on the p-values, and output the results to a file
    #if fold change is 1.0000000000000002, remove those transcripts

    expression_df = expression_df[expression_df["fold-change"] != 1.0000000000000002]
    expression_df["p-value"] = statsmodels.stats.multitest.multipletests(expression_df["p-value"], method="fdr_bh")[1]

    if min_p_value:
        expression_df = expression_df[expression_df["p-value"] < min_p_value[0]]

    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(f"{output_name}.fdr_bh.tsv", sep="\t")
    return expression_df



expression_df = performTTest(expression_df, df_sample_metadata, min_concentration=0)

#read in the expression dataframe from file
expression_df = pd.read_csv(output_name, sep="\t", index_col=0)
expression_df = performBenjaminiHochbergCorrection(expression_df, output_name)
