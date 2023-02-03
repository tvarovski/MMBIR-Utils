import pandas as pd
import os
import cancer_config as cfg

# set environment variables
cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

#expression_data_path = r"C:\Users\twaro\OneDrive\Desktop\development\expression"
#expression_data_path = os.path.realpath(expression_data_path)
expression_data_path_root = cfg.settings["expression_data_path_root"]
expression_data_path = f"{expression_data_path_root}/{username}/TCGA-{cancer}/expression"

#sample_metadata_path = r"C:\Users\twaro\OneDrive\Desktop\development\TCGA-BRCA-RNA-Seq-TSV-manifest.tsv"
metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

#output_name = r"C:\Users\twaro\OneDrive\Desktop\development\expression\expression_data.tsv"
output_name = "expression_data_{cancer}.tsv"

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


def loadSampleMetadata(sample_metadata_path):
    # read in the sample metadata
    sample_metadata = pd.read_csv(sample_metadata_path, sep="\t")

    # create a new column to store the case ID
    sample_metadata["case_id"] = sample_metadata["cases.0.case_id"]

    # return the sample metadata
    return(sample_metadata)


def addCaseIDtoExpressionDataframe(expression_df, expression_metadata):

    # load expression metdata
    expression_metadata=loadSampleMetadata(metadata_location)

    # iterate through the sample_names in expression_df and add the caseID to the expression_df
    for sample_name in expression_df["sample_name"]:
        expression_df.loc[expression_df["sample_name"] == sample_name, "case_id"] = expression_metadata[expression_metadata["file_name"] == sample_name.split("/")[-1]]["case_id"].values[0]


    # save the expression data
    return expression_df


def main():
    # create the expression dataframe
    expression_df = createExpressionDataframe(expression_data_path)

    # add the case ID to the expression dataframe
    expression_df = addCaseIDtoExpressionDataframe(expression_df, loadSampleMetadata(metadata_location))

    # save the expression data
    expression_df.to_csv(output_name, sep="\t", index=False)

if __name__ == "__main__":
    main()



