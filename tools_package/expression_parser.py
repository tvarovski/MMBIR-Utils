import pandas as pd
import os
import cancer_config as cfg

# set environment variables
cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

expression_data_path_root = cfg.settings["expression_data_path_root"]
expression_data_path = f"{expression_data_path_root}/{username}/TCGA-{cancer}/expression"

#metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
#metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

expression_metadata_file = f"TCGA-{cancer}-WXS-expression-metadata.tsv"
expression_metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{expression_metadata_file}"

output_name = f"expression_data_{cancer}.pickle"

# this function is used to extract the expressoion data from the file specified by the path
def processSample(sample_path):

    # read in the sample file
    expression_sample=pd.read_csv(sample_path, sep="\t", comment="#")

    # extract the gene expression data, transpose the data, and remove the first column
    expression_transpose = expression_sample[["gene_name","tpm_unstranded"]].dropna().T
    expression_transpose.columns = expression_transpose.iloc[0]

    # remove the first row (the header) to have gene names as the index
    expression_transpose.drop(expression_transpose.index[0], inplace=True)

    # add the sample name to the dataframe as a new column and reset the index
    expression_transpose["sample_name"] = sample_path
    expression_transpose.reset_index(drop=True, inplace=True)

    # return the dataframe
    return(expression_transpose)


def createExpressionDataframe(expression_data_path):
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


def loadSampleMetadata(expression_metadata_location):

    # read in the sample metadata
    sample_metadata = pd.read_csv(expression_metadata_location, sep="\t")

    # create a new column to store the case ID
    sample_metadata["case_id"] = sample_metadata["cases.0.case_id"]

    # return the sample metadata
    return(sample_metadata)


def addCaseIDtoExpressionDataframe(expression_df, expression_metadata_location):

    # load expression metdata
    expression_metadata = loadSampleMetadata(expression_metadata_location)

    #create sample_name_file column with path removed from sample_name
    expression_df["sample_name_file"] = expression_df["sample_name"].str.split("/").str[-1]
    print("expression_df['sample_name_file']")
    print(expression_df["sample_name_file"].head(10))

    print("expression_metadata['file_name']")
    print(expression_metadata["file_name"].head(10))

    # add the case_id from expression_metadata to the expression_df by matching the file_name
    # if the file_name is not found, the case_id will be NaN
    expression_df["case_id"] = expression_df["sample_name_file"].map(expression_metadata.set_index("file_name")["case_id"])

    # save the expression data
    print(expression_df.head())
    return expression_df


def main():
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

if __name__ == "__main__":
    main()

