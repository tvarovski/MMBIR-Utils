##!/Users/twarowski/mambaforge/bin/python

# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents: 
## filterByDiffExpr
## addMMBToDF
## addMetadataToDF
## perform_UMAP
## plot_UMAP_numerical
## plot_UMAP_categorical

import os
import logging
import umap
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import cancer_config as cfg
from tools import time_elapsed, fancy_status


#parameters#################

cancer = cfg.settings["TCGA-PROJECT"]
pickle_path = f"expression_data_{cancer}.pickle"
df_consolidated_path = f"consolidated_results_{cancer}.tsv"
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

#df_diff_expr_path = f"outputs/ttest_results_{cancer}_minconc0_bh_corrected.tsv"
df_diff_expr_path = f"outputs/ttest_results_BRCA_minconc0_low0.6_high0.2_bh_corrected.tsv"

sample_type = "Primary Tumor" #where to get MMBIR counts from
n_neighbors = 50
min_dist = 0.1
n_components = 2
metric = "canberra"

filter_transcriptome = True
min_pval = 0.05
fold_change_lower = 0.75
fold_change_upper = 1.5

scale_upper = 700
scale_lower = 50

show_plot = True

#############################

logging.basicConfig(level=logging.DEBUG)

def filterByDiffExpr(df, df_diff_expr, min_pval=0.1, fold_change_lower=0.75, fold_change_upper=1.5):

    # take only the genes where the p-value is less than min_pval
    df_diff_expr = df_diff_expr[df_diff_expr["p-value"] < min_pval]
    logging.info(f"There are now {df_diff_expr.shape[0]} differentially expressed transcripts with p-value < {min_pval}.")

    # take only the genes where the fold change is outside fold_change_lower and fold_change_upper
    df_diff_expr = df_diff_expr[(df_diff_expr["fold-change"] < fold_change_lower) | (df_diff_expr["fold-change"] > fold_change_upper)]
    logging.info(f"There are now {df_diff_expr.shape[0]} differentially expressed transcripts with fold-change < {fold_change_lower} or > {fold_change_upper}.")

    #get the gene_id of the differentially expressed genes
    diff_expr_genes = df_diff_expr["gene_id"].tolist()

    #add columns to keep to the list
    diff_expr_genes.extend(['expression_file_name', 'sample_name_file', 'case_id'])

    #drop the columns from the df that are not in the list
    df = df[diff_expr_genes]

    return df

def addMMBToDF(df, df_consolidated, sample_type):

    # filter the df_consolidated dataframe to only include primary tumors
    df_consolidated = df_consolidated[df_consolidated["Sample_Type"] == sample_type]

    #group df_consolidated by Case_ID and get the highest value of Raw_Count and Filtered_Count
    grouping_dict = {"Raw_Count": "max", "Filtered_Count": "max"}
    df_consolidated = df_consolidated.groupby("Case_ID").agg(grouping_dict).reset_index()

    #add a column with Raw_Count to the df by merging on the Case_ID column (in df_consolidated) and the case_id column (in df)
    df = df.merge(df_consolidated[["Case_ID", "Raw_Count", "Filtered_Count"]], left_on="case_id", right_on="Case_ID")
    logging.debug(f"Added the Raw_Count and Filtered_Count columns to the dataframe.")

    #remove the Case_ID column
    df = df.drop(columns=["Case_ID"])

    return df

def addMetadataToDF(df, df_metadata, columns, sample_type=sample_type):

    '''Should drop the duplicates based on filename, not case_id...'''

    #keep only the metadata pertaining to the "Primary Tumor" sample type
    df_metadata = df_metadata[df_metadata["cases.0.samples.0.sample_type"] == sample_type]

    #remove duplicate cases.0.case_id rows, say how many were removed
    logging.info(f"df_metadata shape before removing duplicates: {df_metadata.shape}")
    df_metadata = df_metadata.drop_duplicates(subset=["cases.0.case_id"])
    logging.info(f"df_metadata shape after removing duplicates: {df_metadata.shape}")
    
    #check that there are the same number of rows before and after the merge
    logging.info(f"df shape before merge: {df.shape}")

    df = df.merge(df_metadata[["cases.0.case_id"]+columns], left_on="case_id", right_on="cases.0.case_id")
    logging.debug(f"Added the {columns} columns to the dataframe.")
    df = df.drop(columns=["cases.0.case_id"])

    logging.info(f"df shape after merge: {df.shape}")
    
    return df

@fancy_status
@time_elapsed
def perform_UMAP(rnaseq_data, n_neighbors, min_dist, n_components, metric):

    logging.debug("rna_seq_data shape:")
    logging.debug(rnaseq_data.shape)

    # first, create a umap object, and set the parameters including random_state
    reducer = umap.UMAP(        
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric=metric,
            random_state=42)

    logging.info("Scaling the data")
    scaled_rnaseq_data = StandardScaler().fit_transform(rnaseq_data)
    logging.info("Done scaling the data")

    # then, fit the umap object to the data
    logging.info("Fitting the umap object to the reducer")
    reducer.fit(scaled_rnaseq_data)
    logging.info("Done fitting the umap object to the reducer")

    # then, transform the data using the umap object
    logging.info("Performing UMAP...")
    embedding = reducer.transform(scaled_rnaseq_data)

    logging.debug(f"The shape of the embedding is {embedding.shape}")
    logging.info("Done performing UMAP, returning the embedding")

    return reducer, embedding

def plot_UMAP_numerical(embedding, df, numerical_col, scale_upper=None, scale_lower=None, save_path=None, show_plot=True):

    # plot the umap projection, coloring by the Raw_Count column, set the highest value of color to 1000
    plt.figure(figsize=(12,12))
    # assign marker by figo_stage
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df[numerical_col], s=16, alpha=0.8, cmap='viridis', vmin=scale_lower, vmax=scale_upper)
    plt.title(f'UMAP of RNAseq by {numerical_col}', fontsize=24, fontweight='bold')
    plt.xlabel("UMAP 1", fontsize=18)
    plt.ylabel("UMAP 2", fontsize=18)
    plt.axis('square')
    plt.colorbar(label=f"{numerical_col}")
    #save the figure
    if save_path is not None:
        plt.savefig(f"{save_path}.png", dpi=600)
        logging.info(f"Saved the figure to {save_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_UMAP_categorical(embedding, df, categorical_col, save_path=None, show_plot=True):

    color_labels = df[categorical_col].unique() # get the unique values of the column
    rgb_values = sns.color_palette("tab10", len(color_labels)) # create a list of rgb values
    color_map = dict(zip(color_labels, rgb_values)) # create a dictionary of the color labels and the rgb values
    logging.debug(f"The color labels are {color_labels}")

    plt.figure(figsize=(12,12))
    # assign color by column_var which is a string, so it will be categorical
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df[categorical_col].map(color_map), s=16, alpha=0.5)
    plt.axis('square')

    plt.title(f'UMAP of RNAseq by {categorical_col}', fontsize=24)
    plt.xlabel("UMAP 1", fontsize=18)
    plt.ylabel("UMAP 2", fontsize=18)

    # reuse the color_map for the legend
    patches = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(f"{save_path}", dpi=600)
        logging.info(f"Saved the figure to {save_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_UMAP_custom_IDs(embedding, df, custom_IDs=list, save_path=None, show_plot=True):

    plt.figure(figsize=(12,12))

    #set color to True if the case_id is in the custom_IDs list, otherwise set it to False
    df["color"] = df["case_id"].isin(custom_IDs)
    color_labels = df["color"].unique() # get the unique values of the column
    rgb_values = sns.color_palette("tab10", len(color_labels)) # create a list of rgb values
    color_map = dict(zip(color_labels, rgb_values)) # create a dictionary of the color labels and the rgb values
    logging.info(f"The color map is {color_map}")

    # assign marker by color
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df["color"].map(color_map), s=16, alpha=0.8)
        
    plt.title(f'UMAP of RNAseq by CustomIDs', fontsize=24, fontweight='bold')
    plt.xlabel("UMAP 1", fontsize=18)
    plt.ylabel("UMAP 2", fontsize=18)
    plt.axis('square')
    
    #add a legend
    patches = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    #save the figure
    if save_path is not None:
        plt.savefig(f"{save_path}.png", dpi=600)
        logging.info(f"Saved the figure to {save_path}")
        
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_UMAP_embedding_plotly(embedding, save_path, show=False):

    import plotly.express as px

    fig_2d = px.scatter(
        embedding, x=0, y=1,
        color=df.Filtered_Count,
        color_continuous_scale=px.colors.sequential.Viridis,
        title="UMAP of RNAseq by Filtered_Count",
        hover_name=df.case_id,
        hover_data=[df.Filtered_Count, df.Raw_Count],
        labels={"0": "UMAP 1",
                "1": "UMAP 2",
                "hover_data_0": "Filtered_Count",
                "hover_data_1": "Raw_Count"

                }
    )

    if show:
        fig_2d.show()

    #save the figure to html
    fig_2d.write_html(f"{save_path}")

if __name__ == "__main__":

    # load the data from a pickle file, the consolidated results file, and the differentially expressed genes file
    df = pd.read_pickle(pickle_path)
    df_consolidated = pd.read_csv(df_consolidated_path, sep="\t")
    df_diff_expr = pd.read_csv(df_diff_expr_path, sep="\t")
    df_metadata = pd.read_csv(metadata_location, sep="\t")

    if filter_transcriptome:

        df = filterByDiffExpr(df, df_diff_expr, min_pval=min_pval, fold_change_lower=fold_change_lower, fold_change_upper=fold_change_upper)

    # check if df contains duplicate case_id values, if so need to remove them
    if df["case_id"].duplicated().any():
        logging.info("Duplicate case_id values found, removing duplicates...")
        logging.info(f"The duplicate case_id values are {df[df['case_id'].duplicated()]}")
        df = df.drop_duplicates(subset=["case_id"])
        # check if the duplicates were removed
        if df["case_id"].duplicated().any():
            logging.critical("Duplicate case_id values still found, exiting...")
            exit()
        else:
            logging.info("Duplicate case_id values removed successfully")

    df = addMMBToDF(df, df_consolidated, sample_type)

    # add the specified metadata to the dataframe
    column_var_list_numerical = [   
                                "cases.0.samples.0.days_to_collection",
                                "cases.0.diagnoses.0.age_at_diagnosis",
                                "cases.0.demographic.days_to_birth"
                                ]
    
    column_var_list_categorical = [
                                "cases.0.disease_type",
                                "cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration",
                                "cases.0.diagnoses.0.prior_malignancy",
                                "cases.0.diagnoses.0.prior_treatment",
                                "cases.0.diagnoses.0.tumor_grade",
                                "cases.0.samples.0.preservation_method",
                                "cases.0.samples.0.sample_type", # this is the same as, check both if needed here
                                "cases.0.samples.0.sample_type_id", # this
                                "cases.0.samples.0.submitter_id",
                                "cases.0.submitter_id",
                                "cases.0.diagnoses.0.ajcc_pathologic_stage",
                                "cases.0.diagnoses.0.ajcc_pathologic_m",
                                "cases.0.diagnoses.0.ajcc_pathologic_n",
                                "cases.0.diagnoses.0.ajcc_pathologic_t",
                                "cases.0.diagnoses.0.figo_stage"
                                ]

    column_var_list=column_var_list_numerical + column_var_list_categorical

    try:
        #check if items of column_var_list are present in df_metadata, if not, remove them from column_var_list
        column_var_list_copy = column_var_list.copy()
        for item in column_var_list_copy:
            if item not in df_metadata.columns:
                logging.warning(f"{item} is not in df_metadata, removing it from column_var_list...")
                column_var_list.remove(item)

        df = addMetadataToDF(df, df_metadata, column_var_list)

    except Exception as e:
        logging.error(f"Error adding metadata to dataframe: {e}")
        

    #prepend "Raw_Count" and "Filtered_Count" to column_var_list_numerical
    column_var_list_numerical = ["Raw_Count", "Filtered_Count"] + column_var_list_numerical
    trailing_columns = 5 + len(column_var_list) # add the number of metadata columns to trailing_columns

    #print last trailing_columns + 3 extra columns
    logging.debug(df.iloc[:, -(trailing_columns+3):])
    logging.debug(f"There are {len(df.columns)} columns in the dataframe and {trailing_columns} are trailing")

    # create a list of the columns that we want to use for the umap projection
    #it is all columns except the last few columns which are not part of the data
    rnaseq_data = df.iloc[:, :-trailing_columns]

    # perform the UMAP projection
    reducer, embedding = perform_UMAP(rnaseq_data, n_neighbors, min_dist, n_components, metric)

    # add the Umap axis to the dataframe as columns
    df["UMAP_1"] = embedding[:, 0]
    df["UMAP_2"] = embedding[:, 1]

    # print case_id of points that are less than 5 in the UMAP_1 axis and less than 12 in the UMAP_2 axis
    # save those points to a csv file (just case_id and UMAP_1 and UMAP_2)
    df_UMAP_1_UMAP_2 = df[["case_id", "UMAP_1", "UMAP_2"]]
    df_UMAP_1_UMAP_2 = df_UMAP_1_UMAP_2[(df_UMAP_1_UMAP_2["UMAP_1"] < 5) & (df_UMAP_1_UMAP_2["UMAP_2"] < 12)]
    df_UMAP_1_UMAP_2.to_csv(f"outputs/umap_1_2_{cancer}.csv", index=False)

    #exit()

    # plot the UMAP embedding
    plot_UMAP_embedding_plotly(embedding, show=True, save_path=f"outputs/plots/umap/umap_{cancer}_by_case_id.html")

    # create a folder in the outputs folder to store the plots
    if not os.path.exists("outputs/plots/umap"):
        logging.info("Creating outputs/plots/umap folder")
        os.makedirs("outputs/plots/umap")


    for column_var in column_var_list_numerical:

        if (column_var == "Raw_Count") or (column_var == "Filtered_Count"):

            scale_upper_pass = scale_upper
            scale_lower_pass = scale_lower

        else:
            scale_upper_pass = None
            scale_lower_pass = None

        try:
            plot_UMAP_numerical(embedding, 
                                df, 
                                column_var, 
                                scale_upper=scale_upper_pass, 
                                scale_lower=scale_lower_pass, 
                                save_path=f"outputs/plots/umap/umap_{cancer}_by_{column_var}.png", 
                                show_plot=show_plot)
        except Exception as e:
            logging.error(f"Could not plot numerical_col {column_var}")
            logging.error(e)
            pass


    try:
        # triple:
        custom_IDs=['a5b44d66-c162-46b5-9df2-86305f0385c5', 'a9bb8159-32f0-454c-a946-b3286a52b9d5', '1c3610f7-e0aa-48d7-9a27-0dbaf6e244f9', 'eb2dbb4f-66b6-4525-8323-431970f7a64e', '324bcba2-f6a4-45a6-807c-215bdffcca21', '99e32c59-aa73-43fa-88c9-399ddadb2c72', '18d35983-ea6a-4b70-a209-9bef37595956', '95c53f69-0f05-4348-822a-571f3d757001', '128d198e-9b22-427c-90db-3714455f3a17', '84a88e8d-1e65-4b88-aae6-bf4c7a7e4c33', '295cf595-4a29-46a3-a0c5-e5f08f947031', 'C72CB184-462D-4009-9CDB-848782FF8A76', '53886143-C1C6-40E9-88E6-E4E5E0271FC8', 'DEBA32E4-0E68-4711-941B-3B63BD965AFB', '3afa1e93-1df8-4e4c-aaa4-557463f4bb77', '5ed024e8-d05e-4c65-9441-eda9930ccc82', 'DBF6F981-15CC-40AD-91AC-A66360405FBD', '57AF5C72-0D60-4A6B-B1B4-EC6DAB90F80F', '98E709E7-E195-4B37-9537-F6081AFFB609', 'E8EA576B-EB16-44F7-B241-DAEFC1375388', 'b6633c36-7ce4-4b69-9bf6-30b64d46c66f', 'cf9db1af-17f0-490e-8139-142bd704763a', '6412854a-f874-469e-9d1f-8bd3ae5bd41d', '27b05b15-a44b-45ed-a6e3-e7d1ca488ea9', '6FA2A667-9C36-4526-8A58-1975E863A806', '67C5F371-3FA9-47C5-8B15-C2DD9ACC8519', 'D9FD2724-7DB0-4AF3-AC14-217BDFA5203F', '8183F0FB-2303-4D7B-BCCD-55E5031FC7DF', '2B36853F-34D3-47C5-BA6A-E5A93233D2B1', '49717f75-0f2d-4e1c-9a12-f1cd7877b80a', '3e9f93c0-aa79-4b4c-bd6c-b3325912362a', 'ebc460c2-df88-4dc0-a2c0-aca8072b75ad', '05506f4c-e701-4a9d-ae06-97f066aade43', '96312510-c126-485d-8109-ed81844a1dc3', 'a1da9db1-db11-4f4d-a249-c738db81b87e', '0685edd2-ce1c-4e0e-8dda-393139af4223', 'a9b7d7fe-be31-4f71-afee-c1bfdf511888', 'e9f4f373-37a5-48ad-a1a0-b0d47820111a', 'd5d8e76e-2f2a-49bc-ba0a-d51283611b16', 'b094e8b2-8ece-4b36-8025-18073a8b873c', '016caf42-4e19-4444-ab5d-6cf1e76c4afa', 'c70ee5e1-4703-4996-bb5c-f4cca0fb53ba', '495f9eda-b9be-4eef-bae1-766a4bfbd68e', 'b7f74ae1-6f58-447c-be50-a7666eb19d9a', '9938ce5c-e74e-446e-a932-f096f85cc3b1', '972447aa-4332-47e7-bddd-2eb699dbb664', 'dfe6db17-cf45-488e-bd3c-b8433d7343ca', '17d9e646-6ab3-40b3-a0bc-2c834d3c3213', '71f97b63-c970-44ee-98c2-e02e663d5a40', 'eda6d2d5-4199-4f76-a45b-1d0401b4e54c', '5fd9552a-c742-4388-940d-295d1107ae00', 'd3d545b3-457f-4389-821f-704cb24aff7f', 'D093173F-08AB-4138-BF3C-399C45A6E163', '3AF31FCF-AD0C-4FD9-A8E3-10F9176B5E9D', '9434687A-197C-4959-B6B8-9C05F1DD7F53', 'D9A1C06F-7B50-46B0-878E-89E8C31863AE', '75B3FE55-1A63-426E-867E-2EF52F54778D', '18eb4dfc-556f-4bf3-a411-4780209ed1e0', '5a57dc25-d252-4b22-b192-4de630d7002f', '4f608829-ffc4-4527-886e-6bc764ab29f5', '92b5de82-0221-4df1-8094-80f40c0bb4fa', '20e8106b-1290-4735-abe4-7621e08e3dc8', 'e828455d-4680-41b2-8a1c-1582e3790d62', '786e8dbe-442e-4551-87b3-b4c333b04dd4', '747083ff-0703-431b-aad4-f2adff739516', '6b960b58-28e1-41c6-bd6e-7e669c6aa4ef', '9d166970-07c8-4ca3-9cfa-ed0049df9ecc', 'ad18820b-a804-49c0-ba8a-86c09fa6bce2', 'b97bf89a-7a85-4eef-ae7e-f787aead1f0a', '2DC1CEC9-925A-417F-9E21-3C2143E711B4', '7e1673f8-5758-4963-8804-d5e39f06205b', 'c49e3b18-fd88-48f4-8b01-300692ceb367', 'a855c228-a263-44df-87a0-cbc32187e3f5', '70ab4f23-23c4-409f-a5cc-18a010c3a24e', '4d51159f-ab2f-40e3-a363-847c3654431e', 'f0d8a1fe-e313-44f1-99cc-b965cbeeff0e', '959FF069-8A49-4C9B-85C2-5291CAC0ACFF', '5a17dcd9-5ced-4a69-8069-23c7fd0649d1', '1c40b84e-a0e3-429f-a48c-21566cf881c0', 'fc18d029-9be2-4fa0-9aef-6d647dc55f0b', 'dcd5e079-813a-4c1a-b320-0931468d2bbc', 'aa4244a8-0454-4247-a1c3-357fd51746fa', '97943d87-fed7-4f14-a0a7-c5bfee64c392', '3886b4fb-ba22-4a50-b11a-a6893951f170', '1549dc64-3dab-43fc-96e9-b07d520957e1', '0dca98b0-f43e-45b6-9a02-00092c78678c', 'ae65baeb-6b78-492a-8c63-bb7e93e83dc2', '752ad011-79e0-494f-9868-98bf6feb28f8', '4fba3deb-db94-44b5-a0f2-a575d270779e', '60df7543-6da5-4c75-943b-5800c1e08234', '4ac693e9-10f3-46b3-9d46-df3af7b0d259', 'a6edb6ca-ae9f-4da7-8ebe-92d83d2987fb', '791c5768-f0f5-4ab6-86eb-998e5c4b49e3', 'd8ecba6a-9fff-4993-a799-9a8d2aea524e', '2a84997d-ccee-4f46-bea2-752534f26416', '91a2f2af-e4b1-4a0b-ab20-6a36ce63c533', 'E3935CE4-64D3-4A66-BA11-D308B844B410', 'ac68d219-5670-4ddd-8df6-8aa7ad59e5c7', 'f55dd73d-8c36-440b-84e5-9aae53107775', '521e2140-ae5b-456f-8699-97398d009687', '359f12f9-5c41-48a4-85bc-fd7e307bf7d8', '88db1340-e4bf-451a-87c0-6e9168296f5e', '4da999a0-ef41-4a0b-b1d1-446b39cc855a', 'D9DC3B59-613D-469E-8B4F-6C5A557EB26A', '5C59028F-B8FA-4811-8314-BE3EAED5F364', '9F6BE944-83DE-42AB-8738-F0022F475E61', '6E6B7742-A562-490B-BB72-04A5653852E4', 'D6F7AFC0-1558-43AD-ACB1-2B5311ED2264', '5700C1B3-922A-4401-8AB8-89028908D696', '124B693C-77DC-4FA8-B703-54E8B5054A92', '69FC24FC-CFE2-487A-935F-1A954B30B709', 'B8AEFC48-4A6E-4254-A57F-5F688399B582', 'EA645243-DF49-4466-A255-9F3D4321E357', '398FB71B-CA83-44E7-BF0D-B1CA464B0283', '23C31C2E-336C-4878-A476-CF8D811B4875', 'E783E518-C1E5-4EAC-8EAF-E2D65CCD9692']
        
        logging.info(f"Provided {len(custom_IDs)} custom IDs for UMAP plot...")
        plot_UMAP_custom_IDs(embedding, df, custom_IDs, save_path=f"outputs/plots/umap/umap_{cancer}_by_triple_negative.png", show_plot=True)

    except Exception as e:
        logging.error(f"Error in tripleNeg UMAP plot: {e}")


    for column_var in column_var_list_categorical:

        if column_var == "cases.0.diagnoses.0.figo_stage":
            try:
                #replace the column_var names, specifically remove A, B, C, 1, 2
                df[column_var] = df[column_var].str.replace("A", "")
                df[column_var] = df[column_var].str.replace("B", "")
                df[column_var] = df[column_var].str.replace("C", "")
                df[column_var] = df[column_var].str.replace("1", "")
                df[column_var] = df[column_var].str.replace("2", "")
                logging.info("Removed ('A', 'B', 'C', '1', '2') from figo_stage column values")
            except Exception as e:
                logging.error("Could not remove ('A', 'B', 'C', '1', '2') from figo_stage column values")
                logging.error(e)
                continue
        try:
            plot_UMAP_categorical(embedding,
                                df, 
                                column_var, 
                                save_path=f"outputs/plots/umap/umap_{cancer}_by_{column_var}.png", 
                                show_plot=show_plot)
        except Exception as e:
            logging.error(f"Could not plot categorical_col {column_var}")
            logging.error(e)
            continue