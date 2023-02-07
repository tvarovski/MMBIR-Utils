import pandas as pd
import scipy
import numpy
import statsmodels.stats.multitest
import cancer_config as cfg
import pyensembl
from tools import getCasesAboveMMBThreshold


def performTTest(expression_df, df_sample_metadata, output_name, min_concentration=0):

    #if there are multiple samples associated with a case_id, print a warning
    if len(expression_df["case_id"].unique()) != len(expression_df["case_id"]):
        print("WARNING: There are multiple samples associated with a case_id")

        #print the case_ids that have multiple samples
        df_multiple_samples = expression_df[expression_df.duplicated(subset="case_id", keep=False)]
        print(f"There are {len(df_multiple_samples['case_id'].unique())} cases with multiple samples")

        #ask the user if they want to keep only one sample per case_id
        keep_one_sample = input("Do you want to keep only one sample per case_id? (y/n): ")

        if keep_one_sample == "y":

            print("Keeping only one sample per case_id")
            # keep only one sample per case_id
            expression_df = expression_df.drop_duplicates(subset="case_id", keep="first").copy()

            print(f"There are now {len(expression_df['case_id'].unique())} cases with one sample")
            print(f"Following case_ids kept only the first record: {df_multiple_samples['case_id'].unique()}")
        
        else:

            print("Keeping all samples")
            #keep all samples
            pass
        

    threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_HIGH, min_concentration=min_concentration)
    threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_LOW, below=True, min_concentration=min_concentration)

    # get the IDs of the cases that are high mmbir
    high_mmbir_cases=threshold_mmb_cases_high_df["Case_ID_"].values.tolist()
    low_mmbir_cases=threshold_mmb_cases_low_df["Case_ID_"].values.tolist()


    # Mark the cases that are above the MMBIR threshold as high in the expression dataframe
    expression_df["high_mmbir"] = "none"
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

    #reset index, rename column
    expression_df = expression_df.reset_index()
    expression_df = expression_df.rename(columns={"index": "gene_id"})


    # sort the dataframe by the two-sample t-test, lowest p-value first, and output the results to a file
    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(output_name, sep="\t")
    return expression_df


def performBenjaminiHochbergCorrection(expression_df, output_name, *min_p_value):

    #show column names
    print(expression_df.columns.values.tolist())

    print("Performing Benjamini-Hochberg correction on p-values...")

    #if fold change is 1.0, remove those transcripts
    expression_df = expression_df[expression_df["fold-change"] != 1.0]

    #perform a benjamini-hochberg correction on the p-values, and output the results to a file
    expression_df["p-value"] = statsmodels.stats.multitest.multipletests(expression_df["p-value"], method="fdr_bh")[1]

    if min_p_value:
        expression_df = expression_df[expression_df["p-value"] < min_p_value[0]]

    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(f"{output_name}", sep="\t", index=False)
    return expression_df


def count_calls(func):
    # decorator that counts the number of times a function is called

    #create a wrapper function that counts the number of times the function is called
    def wrapper(*args, **kwargs):
        wrapper.num_calls += 1
        #print out the number of times the function has been called every 1000 times
        if wrapper.num_calls % 1000 == 0:
            print(f"{wrapper.num_calls} calls to {func.__name__}")
        return func(*args, **kwargs)
    
    #initialize the number of calls to 0
    wrapper.num_calls = 0

    return wrapper

@count_calls
def lambdaEnsemblLookup(gene_id, release=104):
    #look up the gene name from the gene ID using pyensembl
    #if the gene ID is not in the database, use the gene ID instead
    try:
        gene_name = pyensembl.EnsemblRelease(release).gene_by_id(gene_id).gene_name
    except Exception as e:
        print(f"couldn't find {gene_id} in the database. Using the gene ID instead.")
        gene_name = gene_id
    return gene_name

def addGeneNameColumnFromGeneID(expression_df, gene_id_column_name):
    #add a gene name column to the expression dataframe by using pyensembl to look up the gene name from the gene ID
    #gene_id_column_name is the name of the column in the expression dataframe that contains the gene ID
    #returns the expression dataframe with the gene name column added

    #apply the pyensembl function to the gene ID column
    try:
        #if the gene ID is in the format ENSG00000157764.13, split it on the period and take the first part
        #if the geneID is not in the database, use the gene ID instead
        expression_df["gene_name"] = expression_df[gene_id_column_name].apply(lambda x: x.split(".")[0])
        expression_df["gene_name"] = expression_df["gene_name"].apply(lambda x: lambdaEnsemblLookup(x))
    except Exception as e:
        print(e)
        print("couldn't find gene name for gene ID. Using gene ID instead")
        expression_df["gene_name"] = expression_df[gene_id_column_name]
    
    return expression_df


if __name__ == "__main__":

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_sample_metadata = pd.read_csv(metadata_location, sep="\t")

    MMBIR_THRESHOLD_LOW = cfg.settings["MMBIR_THRESHOLD_LOW"]
    MMBIR_THRESHOLD_HIGH = cfg.settings["MMBIR_THRESHOLD_HIGH"]

    consolidated_results_path = cfg.settings["consolidated_results_path"]
    expression_df_path = f"expression_data_{cancer}.pickle"
    expression_df = pd.read_pickle(expression_df_path)


    min_concentration=0
    output_name = f"ttest_results_{cancer}_minconc{min_concentration}.tsv"

    expression_df = performTTest(expression_df, df_sample_metadata, output_name, min_concentration=0)

    #read in the expression dataframe from file
    expression_df = pd.read_csv(output_name, sep="\t")
    expression_df = addGeneNameColumnFromGeneID(expression_df, "gene_id")

    output_name = f"ttest_results_{cancer}_minconc{min_concentration}_bh_corrected.tsv"
    expression_df = performBenjaminiHochbergCorrection(expression_df, output_name)
