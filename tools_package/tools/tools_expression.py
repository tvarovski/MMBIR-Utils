import pandas as pd
import scipy
import numpy
import pyensembl
from tools import count_calls, time_elapsed, fancy_status, findThresholdCases
import statsmodels.stats.multitest as smm
import logging

#set the logging level to INFO
logging.basicConfig(level=logging.INFO)

#used by performDiffExprAnalysis
def performExpressionTTest(expression_df, expression_df_metadata, df_sample_metadata, output_name, consolidated_results_path, fraction_high=0.4, fraction_low=0.4, min_concentration=0, investigated_tissue="Primary Tumor"):

    #remove samples that are not Primary Tumor!
    logging.info(f"There are {len(expression_df['case_id'].unique())} cases in the expression_df")
    logging.info(f"Unique sample types: {expression_df_metadata['cases.0.samples.0.sample_type'].unique()}")
    logging.info(f"Keeping '{investigated_tissue}' samples only")
    expression_df_metadata = expression_df_metadata[expression_df_metadata["cases.0.samples.0.sample_type"] == investigated_tissue].copy()
    logging.debug(f"expression_df_metadata: {expression_df_metadata}")
    expression_df = expression_df[expression_df['expression_file_name'].isin(expression_df_metadata["file_name"])].copy()
    logging.info(f"There are now {len(expression_df['case_id'].unique())} cases in the expression_df after removing samples that are not {investigated_tissue}")
    logging.debug(f"expression_df: {expression_df}")

    #if there are multiple samples associated with a case_id, print a warning
    if len(expression_df["case_id"].unique()) != len(expression_df["case_id"]):
        logging.warning("There are multiple samples associated with a case_id")

        #print the case_ids that have multiple samples
        df_multiple_samples = expression_df[expression_df.duplicated(subset="case_id", keep=False)]
        logging.info(f"There are {len(df_multiple_samples['case_id'].unique())} cases with multiple samples")

        #ask the user if they want to keep only one sample per case_id
        keep_one_sample = input("Do you want to keep only one sample per case_id? (y/n): ")

        if keep_one_sample == "y":

            logging.info(f"Keeping only one sample per case_id")
            # keep only one sample per case_id
            expression_df = expression_df.drop_duplicates(subset="case_id", keep="first").copy()

            logging.info(f"There are now {len(expression_df['case_id'].unique())} cases with one sample")
            logging.info(f"Following case_ids kept only the first record: {df_multiple_samples['case_id'].unique()}")
        
        else:
            logging.info("Keeping all samples")
            pass
        
    #old method of getting high+low mmbir cases
    #threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_HIGH, min_concentration=min_concentration)
    #threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, df_sample_metadata, MMBIR_THRESHOLD_LOW, below=True, min_concentration=min_concentration)

    #new method of getting high+low mmbir cases
    df_consolidated = pd.read_csv(consolidated_results_path, sep="\t")

    #default behavior is to filter out cases based on Raw_Counts not Filtered_Counts, but this can be changed by setting filtered=True
    threshold_mmb_cases_high_df, threshold_mmb_cases_low_df = findThresholdCases(df_consolidated, df_sample_metadata, fraction_high=fraction_high, fraction_low=fraction_low, min_concentration=min_concentration, filtered=False)

    #get a list of all the case_ids that are high mmbir

    # get the IDs of the cases that are high mmbir #if using old method, change Case_ID to Case_ID_
    high_mmbir_cases=threshold_mmb_cases_high_df["Case_ID"].values.tolist()
    low_mmbir_cases=threshold_mmb_cases_low_df["Case_ID"].values.tolist()


    # Mark the cases that are above the MMBIR threshold as high in the expression dataframe
    expression_df["high_mmbir"] = "none"
    expression_df.loc[expression_df["case_id"].isin(high_mmbir_cases), "high_mmbir"]="high"
    expression_df.loc[expression_df["case_id"].isin(low_mmbir_cases), "high_mmbir"]="low"

    len_high=len(expression_df[expression_df["high_mmbir"] == "high"])
    len_low=len(expression_df[expression_df["high_mmbir"] == "low"])

    logging.info(f"high cases: {len_high}, low cases: {len_low}")


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
            logging.info(f"transcript_df columns: {transcript_df.columns.values.tolist()}")

            #if there are 2 of the same column name, drop the one whose values are all 0s or NaNs
            if transcript_df.columns.values.tolist().count(transcript) > 1:

                logging.info(f"Dropping one of the columns of {transcript}")
                transcript_df = transcript_df.loc[:,~transcript_df.columns.duplicated()]

                transcript_df[transcript] = transcript_df[transcript].astype(float)
                transcript_df[f"{transcript}_log2"] = transcript_df[transcript].apply(lambda x: numpy.lib.scimath.log2(x+0.1))
            else:
                logging.warning(f"Couldn't drop a column. Skipping {transcript}")
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
        high_mmbir_mean_log2 = high_mmbir_transcript_df[f"{transcript}_log2"].mean()
        low_mmbir_mean_log2 = low_mmbir_transcript_df[f"{transcript}_log2"].mean()

        high_mmbir_mean = high_mmbir_transcript_df[f"{transcript}"].mean()
        low_mmbir_mean = low_mmbir_transcript_df[f"{transcript}"].mean()

        #handle cases where the mean is 0 or extremely low
        if high_mmbir_mean == 0:
            high_mmbir_mean = 0.1
        if low_mmbir_mean == 0:
            low_mmbir_mean = 0.1

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

        output_dict = {"t": t_stat, "p-value": p_value, "fold-change": fold_change}

        logging.debug(f"{transcript}: {output_dict}")

        # store the results and fold change in a dictionary
        expression_dict[transcript] = output_dict


    # create a dataframe from the dictionary where the index is the transcript name and the columns are the results
    expression_df = pd.DataFrame.from_dict(expression_dict, orient="index")

    #remove index, rename column
    expression_df = expression_df.reset_index()
    expression_df = expression_df.rename(columns={"index": "gene_id"})


    # sort the dataframe by the two-sample t-test, lowest p-value first, and output the results to a file
    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(output_name, sep="\t", index=False)

    logging.info(f"Finished writing to {output_name}")

    return expression_df

#used by performDiffExprAnalysis
@count_calls
def lambdaEnsemblLookup(gene_id, release=104):
    #look up the gene name from the gene ID using pyensembl
    #if the gene ID is not in the database, use the gene ID instead
    try:
        gene_name = pyensembl.EnsemblRelease(release).gene_by_id(gene_id).gene_name
    except Exception as e:
        logging.warning(f"Couldn't find {gene_id} in the database. Using the gene ID instead.")
        gene_name = gene_id
    return gene_name

#used by performDiffExprAnalysis
@fancy_status
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
        logging.warning("Couldn't find gene name for gene ID. Using gene ID instead")
        expression_df["gene_name"] = expression_df[gene_id_column_name]
    
    return expression_df

#used by performDiffExprAnalysis
def performBenjaminiHochbergCorrection(expression_df, output_name, *min_p_value):

    logging.debug(f"expression_df columns: {expression_df.columns.values.tolist()}")

    logging.info("Performing Benjamini-Hochberg correction on p-values...")

    #if fold change is 1.0, remove those transcripts
    expression_df = expression_df[expression_df["fold-change"] != 1.0]

    #perform a benjamini-hochberg correction on the p-values, and output the results to a file
    expression_df["p-value"] = smm.multipletests(expression_df["p-value"], method="fdr_bh")[1]

    if min_p_value:
        expression_df = expression_df[expression_df["p-value"] < min_p_value[0]]

    expression_df = expression_df.sort_values(by="p-value")
    expression_df.to_csv(f"{output_name}", sep="\t", index=False)
    return expression_df

#used by performDiffExprAnalysis, can be parallelized
@fancy_status
@time_elapsed
def performDiffExprAnalysis(params):

    cancer = params["cancer"]
    df_metadata = params["df_metadata"]
    #MMBIR_THRESHOLD_LOW = params["MMBIR_THRESHOLD_LOW"]
    #MMBIR_THRESHOLD_HIGH = params["MMBIR_THRESHOLD_HIGH"]
    fraction_high = params["mmbir_fraction_high"]
    fraction_low = params["mmbir_fraction_low"]
    consolidated_results_path = params["consolidated_results_path"]
    expression_df_path = params["expression_df_path"]
    min_concentration = params["min_concentration"]
    outputs_path = params["outputs_path"]
    expression_df_metadata_path = params["expression_df_metadata_path"]
    investigated_tissue = params["investigated_tissue"]

    #read in the expression dataframe from file
    expression_df = pd.read_pickle(expression_df_path)
    expression_df_metadata = pd.read_csv(expression_df_metadata_path, sep="\t")
    output_name = f"{outputs_path}/ttest_results_{cancer}_minconc{min_concentration}_low{fraction_low}_high{fraction_high}.tsv"
    expression_df = performExpressionTTest(expression_df, expression_df_metadata, df_metadata, output_name, consolidated_results_path, fraction_high=fraction_high, fraction_low=fraction_low, min_concentration=min_concentration, investigated_tissue=investigated_tissue)

    #read in the expression dataframe from file
    expression_df = pd.read_csv(output_name, sep="\t")
    expression_df = addGeneNameColumnFromGeneID(expression_df, "gene_id")

    output_name = f"{outputs_path}/ttest_results_{cancer}_minconc{min_concentration}_low{fraction_low}_high{fraction_high}_bh_corrected.tsv"
    expression_df = performBenjaminiHochbergCorrection(expression_df, output_name)

#used by expressionParser
def processSample(sample_path):
    # this function is used to extract the expressoion data from the file specified by the path

    # read in the sample file
    expression_sample=pd.read_csv(sample_path, sep="\t", comment="#")

    # extract the gene expression data, transpose the data, and remove the first column
    expression_transpose = expression_sample[["gene_id", "tpm_unstranded"]].dropna().T
    expression_transpose.columns = expression_transpose.iloc[0]

    # remove the first row (the header) to have gene names as the index
    expression_transpose.drop(expression_transpose.index[0], inplace=True)

    # add the sample name to the dataframe as a new column and reset the index
    expression_transpose["expression_file_name"] = sample_path.split("/")[-1]
    expression_transpose.reset_index(drop=True, inplace=True)

    # return the dataframe
    return(expression_transpose)

#used by expressionParser
@fancy_status
@time_elapsed
def createExpressionDataframe(expression_data_path):
    import os

    # create a new dataframe to store the expression data
    expression_df = pd.DataFrame()

    # find all the sample files in the specified directory
    expression_data = os.listdir(expression_data_path)


    # loop through the sample files, extract the expression data, and append to the dataframe
    for sample in expression_data:
        sample_path = os.path.join(expression_data_path, sample)
        if os.path.isdir(sample_path):
            sample_files=os.listdir(sample_path)
            for sample_file in sample_files:

                # check if the file is a tsv file, and if it is, add it to the dataframe
                if sample_file.endswith(".tsv"):
                    sample_df = processSample(os.path.join(sample_path, sample_file))

                    #concatenate the dataframes
                    expression_df = pd.concat([expression_df, sample_df], ignore_index=True)

    return(expression_df)

#used by expressionParser
def loadSampleMetadata(expression_metadata_location):

    # read in the sample metadata
    sample_metadata = pd.read_csv(expression_metadata_location, sep="\t")

    # create a new column to store the case ID
    sample_metadata["case_id"] = sample_metadata["cases.0.case_id"]

    # return the sample metadata
    return(sample_metadata)

#used by expressionParser
@fancy_status
@time_elapsed
def addCaseIDtoExpressionDataframe(expression_df, expression_metadata_location):

    # load expression metdata
    expression_metadata = loadSampleMetadata(expression_metadata_location)

    #create sample_name_file column with path removed from sample_name
    expression_df["sample_name_file"] = expression_df["expression_file_name"].str.split("/").str[-1]

    logging.debug("expression_df['sample_name_file'] (head)")
    logging.debug(expression_df["sample_name_file"].head(10))

    logging.debug("expression_metadata['file_name'] (head)")
    logging.debug(expression_metadata["file_name"].head(10))

    # add the case_id from expression_metadata to the expression_df by matching the file_name
    # if the file_name is not found, the case_id will be NaN
    expression_df["case_id"] = expression_df["sample_name_file"].map(expression_metadata.set_index("file_name")["case_id"])

    # save the expression data
    logging.debug(expression_df)
    logging.info("Added case_id to expression dataframe")
    
    return expression_df

#used by expressionParser
@fancy_status
@time_elapsed
def expressionParser(params):

    # set environment variables
    expression_data_path = params["expression_data_path"]
    expression_metadata_location = params["expression_metadata_location"]
    output_name = params["output_name"]


    # create the expression dataframe
    expression_df = createExpressionDataframe(expression_data_path)

    # save the expression data, in pickle format
    expression_df.to_pickle(f"{output_name}")

    #read expression df from file
    expression_df = pd.read_pickle(f"{output_name}")

    # add the case ID to the expression dataframe
    expression_df = addCaseIDtoExpressionDataframe(expression_df, expression_metadata_location)

    # overwrite the expression data file
    expression_df.to_pickle(f"{output_name}")